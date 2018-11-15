% h = bemobil_plot_patterns(patterns, chanlocs, varargin)
%
%       Plots patterns, optionally resizing them based on a given weight.
%       Calls EEGLAB's topoplot for each pattern.
%
% In:
%       patterns - channels x patterns matrix of pattern weights
%       chanlocs - EEGLAB channel location struct
%
% Optional:
%       scale - [min, max] scale to use for all subplots
%       fixscale - whether or not to fix the colour scale to the min and 
%                  max of all patterns (1|0, default 1)
%       colorbar - whether or not to add a colorbar to the figure (0|1,
%                  default 1 for fixscale or a given scale, 0 otherwise)
%       weights - vector of weights. patterns will be zoomed in/out
%                 relative to the absolute value of these weights. 
%                 (default: all ones)
%       titles - cell of pattern subplot titles (default: pattern numbers)
%       camzoom - patterns will be zoomed in/out by this much by default,
%                 allowing a tighter fit (default: 1.15)
%       minweight - minimal weight necessary for patterns to be plotted 
%                   (default: 0)
%
% Out:  
%       h - handle of the generated figure
%
%                    Copyright 2018 Laurens R Krol
%                    Team PhyPA, Biological Psychology and Neuroergonomics,
%                    Berlin Institute of Technology

% 2018-06-13 lrk
%   - Adjusted fixscale scale to always be centred around 0
%   - Adjusted subplot titles to reflect negative weights
%   - Added custom scale option
%   - Added colorbar
% 2018-03-19 First version
% 2018 adjustments by Marius Klug

% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <http://www.gnu.org/licenses/>.


function h = bemobil_plot_patterns(patterns, chanlocs, varargin)

% parsing input
ip = inputParser;

addRequired(ip, 'patterns', @isnumeric);
addRequired(ip, 'chanlocs', @isstruct);

addParamValue(ip, 'scale', [], @isnumeric);
addParamValue(ip, 'fixscale', 1, @isnumeric);
addParamValue(ip, 'colorbar', [], @isnumeric);
addParamValue(ip, 'weights', ones(1, size(patterns, 2)), @isnumeric);
addParamValue(ip, 'titles', num2cell(1:size(patterns, 2)), @iscell);
addParamValue(ip, 'camzoom', 1.15, @isnumeric);
addParamValue(ip, 'minweight', 0, @isnumeric);

parse(ip, patterns, chanlocs, varargin{:})

patterns = ip.Results.patterns;
chanlocs = ip.Results.chanlocs;
scale = ip.Results.scale;
fixscale = ip.Results.fixscale;
plotcolorbar = ip.Results.colorbar;
% weights = abs(ip.Results.weights) ./ max(abs(ip.Results.weights)) .* ip.Results.camzoom;
weights = abs(ip.Results.weights);
this_camzoom = ip.Results.camzoom;
weightsign = sign(ip.Results.weights);
titles = ip.Results.titles;
minweight = ip.Results.minweight;

weights(weights<1e-10) = 1e-10;

weights_too_low = weights < minweight;

patternNumbersToPlot = find(~weights_too_low);

weights(weights_too_low)=[];
patterns(:,weights_too_low)=[];

% standardizing weights
weights = weights./max(weights).*this_camzoom;

% setting scale
if fixscale
    scalemin = -max(abs(patterns(:)));
    scalemax = max(abs(patterns(:)));
    samescale = 1;
    if isempty(plotcolorbar), plotcolorbar = 1; end
elseif ~isempty(scale)
    scalemin = scale(1);
    scalemax = scale(2);
    samescale = 1;
    if isempty(plotcolorbar), plotcolorbar = 1; end
else
    samescale = 0;
    plotcolorbar = 0;
end

% opening figure
h = figure;

% creating subplots
ncols = ceil(sqrt(size(patterns, 2)));
nrows = ceil(size(patterns, 2)/ncols);
fprintf('Plotting %d of %d patterns: ', length(weights), length(weights_too_low));
for p = 1:size(patterns, 2)
    fprintf('%d ', patternNumbersToPlot(p));
    subplot(nrows, ncols, p, 'parent', h);
    
    % adding asterisk to title if weight is negative
    if weightsign(p) == -1, titles{p} = [num2str(titles{p}), '*']; end
    title(titles{patternNumbersToPlot(p)});
    
    % calling topoplot
    if samescale, evalc('topoplot(patterns(:,p), chanlocs, ''electrodes'', ''off'', ''maplimits'', [scalemin, scalemax])');
    else,         evalc('topoplot(patterns(:,p), chanlocs,''electrodes'', ''off'')'); end
    
    % adjusting zoom
    camzoom(weights(p));
end
fprintf('\n');

if plotcolorbar
    % adding colorbar
    lastpos = get(subplot(nrows, ncols, ncols * (nrows-1) + 1, 'parent', h), 'Position');
    colorbar('Position', [lastpos(1) - lastpos(3)/4, lastpos(2), lastpos(3)/10, lastpos(4)]);
end
    
end
