function [conn_corrs, conn_weights, optimal_offest ] = tkCCA(X, Y, stride, num_intervals, kernel, kernelpar, numeig, reg_x, reg_y)
% TKCCA Performs temporal kernal CCA on two datasets.
% Input:	- X, Y: data matrices in row format (data as rows)
%           - stride: The interval size of each time lag (default 1)
%           - num_lags: The number of time shifts to the X dataset (must be even)
%			- kernel: kernel type (gauss, gauss-diag, poly, linear)
%			- kernelpar: kernel parameter value
%			- numeig: number of canonical correlation vectors to retrieve 
%			  (default = 1)
%			- reg_x, reg_y: regularization for X and Y datasets

% Output: TBD

    % Note:  shifting the X matrix affects data frame size
    % since there are then empty values on the edges that would affect
    % correlation. This will result in a total data loss of total stride * num_lags * 2 so extra data
    % on the edges of each epoch is necessary. 
    
    if mod(num_intervals, 2)
        error("An even number of time lags must be provided")
    end
    
    %subtract mean
    X = X-repmat(mean(X),size(X,1),1);
    Y = Y-repmat(mean(Y), size(Y,1),1);
        
    %divide variance
    X = X*sqrt(diag(1./diag(X'*X)));
    Y = Y*sqrt(diag(1./diag(Y'*Y)));
    
    % time shift params
    height = size(X, 1);
    width = size(X, 2);
    buffer_size = stride*num_intervals;

    %Create the shift matrix with copies of the original data
    X_padding = zeros(buffer_size, width);
    X_padded = [X; X_padding];
    X_shifted = repmat(X_padded, 1, num_intervals + 1);
    
    %Shift data 
    for interval = 0:num_intervals
        offset = width*interval;
        X_shifted(:, offset+1:offset+width) = circshift(X_shifted(:, offset+1:offset+width),interval*stride,1);
    end

    % remove x padding and equivilant Y data
    X_shifted = X_shifted(buffer_size + 1: height - buffer_size, :);
    Y_shifted = Y(buffer_size + 1: height - buffer_size, :);
    
    [U, V] = km_kcca(X_shifted, Y_shifted, kernel, kernelpar, reg_x, reg_y, numeig);
    
    X_shifts = reshape(X_shifted, [size(X_shifted,1), size(X, 2), num_intervals+1]);
    conn_corrs = zeros(1, num_intervals + 1);
    conn_weights = zeros(width, num_intervals+1);
    
    for i = 1:num_intervals+1
        %Filters (canonical weights) in the input space of each variable
        Wx = X_shifts(:,:,i)' * U;
        Wy = Y_shifted'  * V;
        
        %Correlatoin calculation
        WxX = Wx' * X_shifts(:,:,i)';
        XWx = X_shifts(:,:,i) * Wx;
        YWy = Y_shifted * Wy;
        WyY = Wy * Y_shifted';

        num = WxX * YWy;
        denom = (WxX * XWx) .* (WyY * YWy);
        conn_corrs(i) = num / denom;
        conn_weights(:,i) = Wx;
    end

    [~, optimal_offest] = max(conn_corrs);
end

