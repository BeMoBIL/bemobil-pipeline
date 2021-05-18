function [U, V] = km_kcca(X, Y,kernel,kernelpar,reg_x, reg_y,numeig)
% KM_KCCA performs kernel canonical correlation analysis.
% Input:	- X, Y: data matrices containing one datum per row
%			- kernel: kernel type (gauss, gauss-diag, poly, linear)
%			- kernelpar: kernel parameter value
%			- reg: regularization
%			- numeig: number of canonical correlation vectors to retrieve 
%			  (default = 1)
% Output:	- U, V:canonical weights for the overall dataset
% USAGE: [y1,y2] = km_kcca(X1,X2,kernel,kernelpar,reg_x, reg_y, numeig)
%
% Author: Steven Van Vaerenbergh (steven *at* gtas.dicom.unican.es), 2012.
%
% The algorithm in this file is based on the following publications:
% D. R. Hardoon, S. Szedmak and J. Shawe-Taylor, "Canonical Correlation
% Analysis: An Overview with Application to Learning Methods", Neural
% Computation, Volume 16 (12), Pages 2639--2664, 2004.
% F. R. Bach, M. I. Jordan, "Kernel Independent Component Analysis", Journal
% of Machine Learning Research, 3, 1-48, 2002.
%
% This file is part of the Kernel Methods Toolbox for MATLAB.
% https://github.com/steven2358/kmbox

N = size(X,1);	% number of data


% handle the number of eigenvectors to return
if nargin<6,
    numeig = 1;
end
if ischar(numeig),
    if strcmp(numeig,'all'),
        numeig = size(X,1); % upper bound, to be reduced later
    end
end

I = eye(N); Z = zeros(N);
N0 = eye(N)-1/N*ones(N);

% get kernel matrices
K1 = N0*km_kernel(X,X,kernel,kernelpar)*N0;
K2 = N0*km_kernel(Y,Y,kernel,kernelpar)*N0;

% restrict number of eigenvalues with rank information
minrank = min(rank(K1),rank(K2));
numeig = min(numeig,minrank);

% 3 GEV options, all of them are fairly equivalent

% % option 1: standard Hardoon
% R = [Z K1*K2; K2*K1 Z];
% D = 1/2*[K1*(K1+reg*I) Z; Z K2*(K2+reg*I)];
% % R = R/2+R'/2;   % avoid numerical problems
% % D = D/2+D'/2;   % avoid numerical problems

% option 2: simplified Hardoon
% R = [Z K2; K1 Z];
% D = [K1+reg*I Z; Z K2+reg*I];
% % R = R/2+R'/2;   % avoid numerical problems
% % D = D/2+D'/2;   % avoid numerical problems

% % option 3: Kettenring-like generalizable formulation
R = 1/2*[K1 K2; K1 K2];
D = [K1+reg_x*I Z; Z K2+reg_y*I];

% solve generalized eigenvalue problem
[alphas,betas] = eig(R,D);
[betass,ind] = sort(real(diag(betas)));
alpha = alphas(:,ind(end:-1:end-numeig+1));
alpha_norms = sqrt((sum(alpha.^2,1)));
alpha = alpha./repmat(alpha_norms,2*N,1);
conn_corrs = betass(end:-1:end-numeig+1);

% expansion coefficients
U = alpha(1:N,:);
V = alpha(N+1:end,:);

end
