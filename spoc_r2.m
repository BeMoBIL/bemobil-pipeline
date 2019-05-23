function [W, A, Cxx, Cxxz, Cxxe, r2_values, r_values, lambda_values, fv] = spoc_r2(X, z, varargin)
% Source Power Co-modulation Analysis (SPoC).
%
% Optimizes spatial filters such that the power of the filtered
% signal maximally correlates with the univariate target function z, as
% described in Dahne et al., 2014a.
%  
% [W, A, Cxx, Cxxz, Cxxe, r2_values, r_values, lambda_values] = ...
%                     spoc_r2(X, z,  <'keyword1', value1, 'keyword2', value2, ...>);
%
% Input:
% X     - epoched data, with size(X) = [n_channels, n_samples_per_epoch_x , n_epochs]
% z     - target function, with size(z) is [n_samples_per_epoch_z, n_epochs]
% --> the number of samples per epoch can differ between X and z
%
% Keyword input:
% 'n_spoc_components' - the number of spoc components to be optimized,
%                           Default is 1. This version of SPOC_r2 will call
%                           itself recursively to generate up to n = rank
%                           components
%    
% 'use_log' - if 1, the correlation of log-power to z is optimized. Default is 0.
% 'Cxxe' - trialwise covariance matrices. dim(Cxxe) = [n_channels, n_channels, n_epochs]
%           If Cxxe is given, then X is ignored! This can make
%           crossvalidation much faster because Cxxe doesn't have to be
%           computed in each fold
%
% More optional keyword input:
% 'n_repeats' - Number of restarts of the optimization, each time with
%                   random initialization. Note that each component gets
%                   this many restarts! Default is 10.
% 'pca_X_var_explained' - Optional dimensioality reduction with PCA, using
%                       the number of PCA components is determined by how
%                       much cumulative variance they should explain. This
%                       number must be between 0 and 1. Default is 1, i.e.
%                       no dimensionality reduction. 
% 'verbose' - If this value is larger than 0, some info will be printed
%               during bootstrapping. Default value is 1.
%
% For more options see the first few lines of the code.... 
% 
% Output:
%   W       - filter matrix, with filters in the columns
%   A       - pattern matrix, with corresponding patterns in the columns
%   Cxx     - data covariance matrix
%   Cxxz    - data covariance matrix weighted by z
%   Cxxe    - trialwise covariance matrix computed from X
%   r2_values   - squared r_values
%   r_values    - vector with r(i) = correlation between the power of the
%                   data after projection on W(:,i) and the target function
%   lambda_values  - vector with lambda(i) = covariance between the power of the
%                   data after projection on W(:,i) and the target function
%
% Example function calls:
%
% The most simple way, using default values for all parameters:
% >> W = spoc(X,z); 
%
% Specifying some parameter values. Here, optimize three components and
% optimize the correlation between log-power of the components and z:
% >> W = spoc(X,z, 'n_spoc_components', 3, 'use_log', 1); 
%
%
% IMPORTANT NOTES:
% - Make sure the 'utils' folder is on the path, because this function
% requires the minFunc optimization package and some helper functions from
% the utils folder
% - Dimensionality reduction with SSD is advised. See Haufe et al. 2014
%
% References:
%
% S. Dahne, F. C. Meinecke, S. Haufe, J. Hohne, M. Tangermann, K. R. Muller, V. V. Nikulin, 
% "SPoC: a novel framework for relating the amplitude of neuronal
% oscillations to behaviorally relevant parameters",
% NeuroImage, 86(0):111-122, 2014
%
% S. Haufe, S. Dahne, V. V. Nikulin, 
% "Dimensionality reduction for the analysis of brain oscillations"
% NeuroImage 101:583-597 2014
%
% sven.daehne@tu-berlin.de, 2014

opt= propertylist2struct(varargin{:});
opt= set_defaults(opt, ...
    'n_spoc_components', 1, ...
    'use_log', 0, ...
    'n_repeats', 10, ...
    'maxIter', 200, ...
    'pca_X_var_explained', 1, ...
    'verbose',1, ...
    'Cxx', [], ...
    'Cxxz', [], ...
    'Cxxe', [],...
    'is_white', 0);

%If X fed in as matrix, transform to cell array of epochs
%Transform matrix dimensions
if not(iscell(X)) && not(isempty(X))
    [~,~, n_eps] = size(X);
    X_mat = X;
    clear X;
  
    X = cell(1,n_eps);
    for epoch = 1:n_eps
        X{epoch} = X_mat(:,:,epoch).'; 
    end
else
    [~,n_eps] = size(X);
    for epoch = 1:n_eps
        X{epoch} = X{epoch}.';
    end
end

if not(isempty(X))
    [~, n_channels] = size(X{1});
    [~, N_e] = size(X);
else
    [~, n_channels,N_e] = size(opt.Cxxe); %needs to be updated
end
if not(size(z,2) == N_e)
    error('X and z must have the same number of epochs!!!');
end

% normalize the target funtion to have zero mean and unit-variance and
% average over all epochs
z = mean(z,1); % take the average within epochs
z = (z-mean(z(:)))./std(z(:));


if opt.is_white
    Cxx = eye(n_channels);
elseif isempty(opt.Cxx)
    % ordinary covariance matrix
    Cxx = cov(vertcat(X{:}));
else
    Cxx = opt.Cxx;
end

if isempty(opt.Cxxe)
    % create trial-wise covariance matrices
    Cxxe = zeros(n_channels, n_channels, N_e);
    for e=1:N_e
        X_e = squeeze(X{e});
        C_tmp = cov(X_e);
        Cxxe(:,:,e) = C_tmp;
    end
else
    Cxxe = opt.Cxxe;
end
if isempty(opt.Cxxz)
    % create the z-weighted covariance matrix
    Cxxz = create_Cxxz(Cxxe, z);
else
    Cxxz = opt.Cxxz;
end

% whiten the data
[M, ~] = whiten_data(opt.pca_X_var_explained, Cxx);
n_components = size(M,1);

Cxxe_white = zeros(n_components, n_components, N_e);
for e=1:N_e
    Cxxe_white(:,:,e) = M * Cxxe(:,:,e) * M';
end
Cxxz_white = M * Cxxz * M';
Cxx_white = M * Cxx * M';

% do SPoC_lambda, i.e. covariance approximation, to get the first starting
% value of w for further numerical optimization
[W, D] = eig(Cxxz_white);
[~, sorted_idx] = sort(abs(diag(D)), 'descend');
W = W(:, sorted_idx);

% optimize the first filter, in whitened space
[w, r2_value] = optimize_W(z, Cxxe_white, W(:,1), opt);


W = w; % start the filter set
r2_values = zeros(1,opt.n_spoc_components);
r2_values(1) = r2_value;
% compute the next components/filters recursively
if opt.n_spoc_components > 1
    opt.n_spoc_components = opt.n_spoc_components - 1;
    opt.pca_X_var_explained = 1; % make sure there is no more dim reduction
    opt.is_white = 1;
    % create a basis of the space that is orthogonal
    B = null(W');
    % project the data to the null-space and thereby reduce the dimensionality
    Cxxe_redux = zeros(n_components-1, n_components-1, N_e);
    for k=1:N_e
        Cxxe_redux(:,:,k) = B' * Cxxe_white(:,:,k) * B;
    end
    opt.Cxxe = Cxxe_redux;
    opt.Cxx = B' * Cxx_white * B;
    opt.Cxxz = B' * Cxxz_white * B;
    % repeat SPoC on the subspace data
    [W_tmp, ~,~,~,~,r2_values_tmp] = spoc_r2([], z, opt);
    % project the filters back to the space of original dimensionality
    W = [W B * W_tmp];
    r2_values(2:end) = r2_values_tmp;
end

W = M'*W; % project back to original (un-whitened) channel space

% Normalize the filters such that the extracted components have unit
% variance.
for k=1:size(W,2)
    W(:,k) = W(:,k) / sqrt(squeeze(W(:,k)'*Cxx*W(:,k)));
end

if nargout > 5
    fv = get_var_features(W, Cxxe);
    R = corrcoef([z', fv']);
    r_values = R(1,2:end);
    C = cov([z', fv']);
    lambda_values = C(1,2:end);
end

A = Cxx * W;



function [W, r2_values] = optimize_W(z, Cxxe, w_spoc, opt)
N_x = size(Cxxe,1);
minFunc_opt = struct('DerivativeCheck', 'off', 'MaxIter', opt.maxIter, ...
    'Display', 'off', 'useComplex', 0, 'numDiff', 0, 'TolFun', 1e-02, 'TolX', 1e-05);
w_tmp = zeros(N_x, opt.n_repeats);
fval_tmp = zeros(1, opt.n_repeats);
warning off
for k=1:opt.n_repeats
    if k==1
        w_start = w_spoc;
    else
        w_start = randn(N_x, 1);
    end
    %[w_tmp(:,k) fval_tmp(k), exitflag, output] = minFunc(@fval_grad, w_start, minFunc_opt, z, Cxxe);
    [w_tmp(:,k), fval_tmp(k)] = minFunc(@fval_grad, w_start, minFunc_opt, z, Cxxe, opt);
end
warning on
[min_r2, idx] = min(fval_tmp);
W = w_tmp(:,idx);
r2_values = abs(min_r2);


function [fval, grad] = fval_grad(w, z, Cxxe, opt)
% Computes the objective function value (squared correlation between z and z_est) and the
% gradient at the current w
[N_x, ~, N_e] = size(Cxxe);
grad_px_wrt_wx = zeros(N_x, N_e);

for e=1:N_e
    grad_px_wrt_wx(:,e) = 2 * Cxxe(:,:,e) * w;
end
px = 0.5* w' * grad_px_wrt_wx;

if opt.use_log
    phi_x = log(px);
    grad_phi_x_wrt_wx = (1./repmat(px, N_x, 1)) .* grad_px_wrt_wx;
else
    phi_x = px;
    grad_phi_x_wrt_wx = grad_px_wrt_wx;
end

mean_phi_x = mean(phi_x);
cov_xz = mean((phi_x-mean_phi_x) .* (z-mean(z)));
sigma_x = std(phi_x);
sigma_z = std(z);
% objective function value -> squared correlation
fval = -1* (cov_xz/(sigma_x*sigma_z))^2;

% gradient
grad_x = gradient_corr_phi_z_wrt_w(phi_x, z, grad_phi_x_wrt_wx);
grad = 2 * fval * grad_x;
     


function Cxxz = create_Cxxz(Cxxe, z)
% create the z-weighted covariance matrix
N_e = size(Cxxe,3);
Z = zeros(size(Cxxe));
for e=1:N_e
    Z(:,:,e) = z(e);
end
Cxxz = sum(Cxxe.*Z , 3) ./ N_e; 

function fv = get_var_features(W, Cxxe)
Ne = size(Cxxe,3);
fv = zeros(size(W,2),Ne);
for e=1:Ne
    fv(:,e) = diag(W'*Cxxe(:,:,e)*W);
end


function grad = gradient_corr_phi_z_wrt_w(phi, z, grad_phi_wrt_w)
%
% phi -     size(phi) = [1, N]
% z -       size(z) = [1, N]
% grad_phi_wrt_w - size(grad_phi_wrt_w) = [d, N]

[d, N] = size(grad_phi_wrt_w);

phi = phi - mean(phi);
z = z - mean(z);
% grad_phi_wrt_w = grad_phi_wrt_w - repmat(mean(grad_phi_wrt_w,2), 1, N);

C = cov(phi, z);
cov_phi_z = C(1,2);
var_phi = C(1,1);
var_z = C(2,2);

delta_z_phi = z - phi*cov_phi_z/var_phi;
grad = mean(repmat(delta_z_phi, d, 1) .* grad_phi_wrt_w, 2) / sqrt(var_phi*var_z);

