function [W, A, lambda, C_s, X_ssd] = ssd(X, freq, sampling_freq, filter_order, epoch_indices)
% SSD - Spatio-Spectral Decomposition, based on Haufe implementation
%
% [W, A, lambda, C_s, X_ssd] = ssd(X, freq, sampling_freq, filter_order, epoch_indices)
%
% This is a function for the extraction of neuronal oscillations 
% with optimized signal-to-noise ratio. The algorithm maximizes 
% the power at the center frequency (signal of interest) while simultaneously suppressing it
% at the flanking frequency bins (considered noise). 
% 
% INPUT: 
%     X -     a matrix of size TxM with T samples and M channels of 
%             raw EEG/MEG/LFP data
%     freq -  3 x 2 matrix with the cut-off frequencies. 
%             First row: cut-off frequencies for band-pass of the to be extracted 
%             oscillations.
%             Second row: cut-off frequencies for the lowest and highest 
%             frequencies defining flanking intervals.
%             Third row: cut-off frequencies for the band-stop filtering of 
%             the central frequency process.
%     sampling_freq -     the sampling frequency (in Hz) of X (the data)
%     filter_order  -     filter order used for butterworth bandpass and 
%                         bandstop filtering. If unsure about it, use [] then
%                         the default value of order = 2 will be used
%     epoch_indices -     a matrix of size N x 2, where N is the number of 
%                         good (i.e. artifact free) data segments. Each row of 
%                         this matrix contains the first and the last sample index of
%                         a data segment, i.e. epoch_indices(n,:) = [1000, 5000]
%                         means that the n'th segment starts at sample 1000 and
%                         ends at sample 5000. If all data is useable or if unsure,
%                         use [], then the default of [1, size(X,1)] will be used. 
%                         
% OUTPUT:
%     W -     the de-mixing matrix. Each column is a spatial filter and the
%             timecourse of the SSD components is extracted with X * W
%     A -     the spatial patterns (also called mixing matrix) with the i'th column
%             corrresponding to the i'th SSD component
%     lambda - the eigenvalues corresponding to each component. The stronger
%              the eigenvalue the better is the ratio between the signal and noise. 
%              The components are sorted in the descending order (first components 
%              have the largest SNR)
%      C_s -  the covariance matrix of X after bandpass filtering with the band
%             defined in freq(1,:)
%      X_ssd - the bandpass filtered data projected onto the SSD components, 
%              i.e. X_ssd = X_s * W, where X_s is the bandpass filtered version of X

% Copyright (c) [2014] [Stefan Haufe, Sven Daehne, Vadim Nikulin]
% 
% Permission is hereby granted, free of charge, to any person obtaining a copy
% of this software and associated documentation files (the "Software"), to deal
% in the Software without restriction, including without limitation the rights
% to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
% copies of the Software, and to permit persons to whom the Software is
% furnished to do so, subject to the following conditions:
% 
% The above copyright notice and this permission notice shall be included in all
% copies or substantial portions of the Software.
% 
% THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
% IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
% FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
% AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
% LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
% OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
% SOFTWARE.


%% check input arguments

%Transpose since EEGLab data struct is stored horizontally
X = X.';

% make sure FREQS has the correct dimensions
if not( size(freq,1)==3 && size(freq,2)==2 )
  error('freq must be a 3 by 2 matrix, i.e. three bands must be specified!');
end

% check the given frequency bands
signal_band = freq(1,:); % signal bandpass band
noise_bp_band = freq(2,:); % noise bandpass band
noise_bs_band = freq(3,:); % noise bandstop band
if not( noise_bs_band(1) < signal_band(1) && ...
        noise_bp_band(1) < noise_bs_band(1) && ...
        signal_band(2) < noise_bs_band(2) && ...
        noise_bs_band(2) < noise_bp_band(2) )
  error('Wrongly specified frequency bands!\nThe first band (signal band-pass) must be within the third band (noise band-stop) and the third within the second (noise band-pass)!');
end

% default values for optional arguments
if isempty(filter_order)
    filter_order = 2;
end
if isempty(epoch_indices)
    epoch_indices = [1, size(X,1)];
end

% indices of good segments
ind=[];
for n=1:size(epoch_indices,1)
    ind=[ind epoch_indices(n,1):epoch_indices(n,2)];
end

%% filtering of data

% Creating filters
[b,a]=butter(filter_order, signal_band/(sampling_freq/2));
[b_f,a_f]=butter(filter_order, noise_bp_band/(sampling_freq/2));
[b_s,a_s]=butter(filter_order, noise_bs_band/(sampling_freq/2),'stop');


% Covariance matrix for the center frequencies (signal)
X_s = filtfilt(b,a,X);
C_s = cov(X_s(ind,:),1);

% Covariance matrix for the flanking frequencies (noise)
X_tmp = filtfilt(b_f,a_f,X);
X_tmp = filtfilt(b_s,a_s,X_tmp);
C_n = cov(X_tmp(ind,:),1);
clear X_tmp


%% Generalized eigenvalue decomposition

% dim-reduction of X does not have full rank
C = C_s;
[V, D] = eig(C);
[ev_sorted, sort_idx] = sort(diag(D), 'descend');
V = V(:,sort_idx);
% compute an estimate of the rank of the data
tol = ev_sorted(1) * 10^-6;
r = sum(ev_sorted > tol);
if r < size(X,2)
    fprintf('SSD: Input data does not have full rank. Only %d components can be computed.\n',r);
    M = V(:,1:r) * diag(ev_sorted(1:r).^-0.5);
else
    M = eye(size(X,2));
end


C_s_r = M' * C_s * M;
C_n_r = M' * C_n * M;
[W,D]= eig(C_s_r,C_s_r+C_n_r);
[lambda, sort_idx] = sort(diag(D), 'descend');
W = W(:,sort_idx);

W = M * W;
% A is the matrix with the patterns (in columns)
A = C * W / (W'* C * W);


%% apply SSD filters to the data

if nargout > 4
    X_ssd = X_s * W;
end
