% bemobil_plot_trail()  - Plots the trail map of the selected eye tracking
%                         data in the EEG data strucutre. Takes both
%                         continuous and epoched data. The continuous data
%                         is not recommended as it may take a long time to
%                         display. The thickness of the trail map
%                         represents the velocity of the movement (thin =
%                         slow, thick = fast).
%
% Usage:
%   >> [trailmap, handleList] = bemobil_plot_trail(EEG, trailFactor, trailColor, trailSize,...
%                               EyePosXName, EyePosYName, norm, figPos)
% Inputs:
%   EEG                     - current EEGLAB EEG structure. Takes only a
%                             single EEG struct at the time.
%   trailFactor             - factor to determine the smoothness of trail (OPTIONAL ARGUMENT)
%   trailColor              - colormap or N x 3 matrix of colors, where N > 2.
%                             Takes also standardized maps, eg. 'summer', 'magma', or
%                             'cividis' (DEFAULT = 'magma') (OPTIONAL ARGUMENT)
%   trailSize               - set the size for the thickness of trail (e.g. [1 30], where [min max])
%                             plot. The thickness represents the velocity
%                             of the movement (DEFAULT = [1 20]) (OPTIONAL ARGUMENT)
%   EyePosXName             - the name of the channel containing the
%                             x-coordinate of the eye-tracking, eg. 'left_x_data'
%   EyePosYName             - the name of the channel containing the
%                             y-coordinate of the eye-tracking, eg. 'left_y_data'
%   norm                    - normalize the data so that it fits between [-1 1] (DEFAULT = 1)
%   figpos                  - set the figure size and position (OPTIONAL ARGUMENT)
%
% Outputs:
%   [trailmap, handleList]  - eye tracking data structure with fieldnames,
%                             and the figure handles when the figures are open.
%
% Example:
%   [trailmap, handleList] = bemobil_plot_trail({ALLEEG(1) ALLEEG(2)},
%                            'EyePosXName', 'CombinedEyes_X',
%                            'EyePosYName', 'CombinedEyes_Y',
%                            'trailColor', 'viridis',
%                            'trailsize', [1 50], 'norm', 1);
%
% See also:
%   EEGLAB, bemobil_plot_eyetracking, bemobil_plot_heatmap, bemobil_plot_motion
%
% Author: Zak Djebbara, April 2022

function [trailmap, handleList] = bemobil_plot_trail(EEG, varargin)

if nargin == 0; help bemobil_plot_trail; return; end

p = inputParser;
epoched = [];
handleList = [];

addRequired(p, 'EEG');
addParameter(p, 'trailFactor', 1, @isnumeric);
addParameter(p, 'trailColor', 'magma');
addParameter(p, 'trailSize', [1 20], @isnumeric);
addParameter(p, 'EyePosXName', 'Eyetracking_Screen_X_both', @ischar);
addParameter(p, 'EyePosYName', 'Eyetracking_Screen_Y_both', @ischar);
addParameter(p, 'figPos', [300 300 600 300], @isnumeric);
addParameter(p, 'norm', 1, @isnumeric);

parse(p, EEG, varargin{:});

trailFactor = p.Results.trailFactor;
trailColor = p.Results.trailColor;
trailSize = p.Results.trailSize;
figPos = p.Results.figPos;
EyePosXName = p.Results.EyePosXName;
EyePosYName = p.Results.EyePosYName;
norm = p.Results.norm;

tiledlayout(figure, 'flow', 'TileSpacing','compact' , 'Padding', 'compact');
for i = 1:length(EEG)
    % Checking if data is epoched or not
    if isempty(EEG{i}.epoch)
        warning(['No epochs detected. Generates continuous plot.']);
        epoched = 0; titleplot = "continuous data";
    elseif ~isempty(EEG{i}.epoch)
        epoched = 1; titleplot = "averaged epoched data";
    end
    pause(.1)

    % Identifying the channel names
    xchanidx = [];
    ychanidx = [];
    for j = 1:length(EEG{i}.chanlocs)
        if strcmp(EEG{i}.chanlocs(j).labels, EyePosXName); xchanidx = [xchanidx, j];end
        if strcmp(EEG{i}.chanlocs(j).labels, EyePosYName); ychanidx = [ychanidx, j];end
    end
    if size(ychanidx,2)> 1 || size(xchanidx,2)> 1; error('More than one channel share the same name for eyetracking channel.');end
    if isempty(ychanidx) || isempty(xchanidx); error('Eyetracking channel could not be found.');end

    % Identifying data
    if epoched 
        xdata       = mean(EEG{i}.data,3);
        xdata       = xdata(xchanidx,:);
        xdatadfdt   = normalize(gradient(xdata),'range',[min(xdata) max(xdata)]);
        ydata       = mean(EEG{i}.data,3);
        ydata       = ydata(ychanidx,:);
        ydatadfdt   = normalize(gradient(ydata),'range',[min(ydata) max(ydata)]);
        tickstep    = length(EEG{i}.times)/5;
        ticksX      = 1:round(tickstep):length(xdata);
        ticksX      = [ticksX, length(xdata)];
        
        dataduration = (round(EEG{i}.times(end))-round(EEG{i}.times(1)));
        durationunit = "ms.";
    elseif ~epoched
        xdata = EEG{i}.data(xchanidx,:);
        ydata = EEG{i}.data(ychanidx,:);
        tickstep = length(EEG{i}.times)/5;
        ticksX = 1:round(tickstep):length(xdata);
        ticksX = [ticksX, length(xdata)];
        
        dataduration = (round(EEG{i}.times(end))-round(EEG{i}.times(1)))/1000;
        durationunit = "s.";
    end
    if norm; xdata = normalize(xdata, 'range', [-1 1]); xdata = normalize(xdata, 'center');end
    if norm; ydata = normalize(ydata, 'range', [-1 1]); ydata = normalize(ydata, 'center');end
    x = interp(xdata,trailFactor);
    y = interp(ydata,trailFactor);
    dfdt = normalize(interp(mean([xdatadfdt;ydatadfdt]),trailFactor), 'range', trailSize);
    
    % Setting up the colormap
    try
        if ischar(trailColor)
            color = eval([trailColor '(' num2str(length(x)) ');']);
        elseif isnumeric(trailColor)
            color = bemobil_makecmap(trailColor, length(x));
        end
    catch
        error(['Colormap not recognized. Please use standard maps, eg. "magma", "summer", or "winter"'...
            newline 'or matrix of at least two rows, eg. [0 0 0; .5 .5 .5];']);
    end
    
    % Plotting 
    nexttile(i);
    scatter(x,y,dfdt,color,'filled');
    colormap(color);
    colorbar(gca, 'Location', 'southoutside','XTickLabel',{'Start','Midway','End'},'XTick',[0 .5 1]);
    title("Trailmap of " + titleplot + "." + newline + "Duration: ~" + dataduration + durationunit);
    xlabel('x coordinate'); ylabel('y coordinate');
    if norm; xlim([-1 1]); ylim([-1 1]);end
    if ~norm; xlim([min(xdata) max(xdata)]); ylim([min(ydata) max(ydata)]);end
    set(gcf, 'Position', figPos);
    axis('image');
    set(gcf, 'Color', 'w');
    handleList = [handleList; gcf];
    
    trailmap{i}.x = x;
    trailmap{i}.y = y;
end
end