% bemobil_makecmap()    - Creates a linearly related colormap based on a (N x 3) matrix.
%                         
%
% Usage:
%   >> [cmap] =   bemobil_makecmap( colormatrix , numcolors)

% Inputs:
%   colormatrix             - [N x 3] color matrix, e.g. 
%                             [255 255 255; 150 150 150; 100 100 100]
%   numcolors               - number of colors to produce.
%
% Outputs:
%   cmap                    - interpolated color matrix
%
% See also:
%   bemobil_plot_eyetracking, bemobil_plot_heatmap, bemobil_plot_patterns
%
% Author: Zak Djebbara, 2022
function [cmap] = bemobil_makecmap(input, numcmap)
if nargin == 0; help bemobil_makecmap; end
if size(input,2) < 3 || size(input,1) < 2; error('Input does not match criteria. Please make sure the input is at least 2 x 3 matrix'); end
if ~isnumeric(numcmap); error('Please make sure the numbcmap variable is a number'); end

cmap = [];

totalnum = numcmap-size(input, 1);
if totalnum<0; error('There are more colors inserted than wished to generate.'); end

tweennum = ceil(totalnum/(size(input,1)-1));
amountoftimes = ceil(totalnum/tweennum);

for i1 = 1:amountoftimes
    newr = []; newg = []; newb = [];
    if size(cmap,1)+tweennum < numcmap
        newr = interp1(input(i1:i1+1,1)', linspace(1,2, tweennum+2));
        newg = interp1(input(i1:i1+1,2)', linspace(1,2, tweennum+2));
        newb = interp1(input(i1:i1+1,3)', linspace(1,2, tweennum+2));
    else
        tweennum = numcmap-size(cmap,1)-1;
        newr = interp1(input(i1:i1+1,1)', linspace(1,2, tweennum+2));
        newg = interp1(input(i1:i1+1,2)', linspace(1,2, tweennum+2));
        newb = interp1(input(i1:i1+1,3)', linspace(1,2, tweennum+2));
    end
    cmap = [cmap; newr', newg', newb'];
    cmap = unique(cmap, 'rows', 'stable');
end
end