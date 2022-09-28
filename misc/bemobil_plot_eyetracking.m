% bemobil_plot_eyetracking()  - Plots the trail map of the selected eye tracking
%                             data in the EEG data strucutre. Takes both
%                             continuous and epoched data. The continuous data
%                             is not recommended as it may take a long time to
%                             display. The thickness of the trail map
%                             represents the velocity of the movement (thin =
%                             slow, thick = fast).
%
% Usage:
%   >> [trailmap, handleList] = bemobil_plot_eyetracking(EEG, trailColor,...
%                               EyePosXName, EyePosYName, names,...
%                               norm, figPos)
% Inputs:
%   EEG                     - current EEGLAB EEG structure. Takes only a
%                             single EEG struct at the time.
%   EyePosXName             - the name of the channel containing the
%                             x-coordinate of the eye-tracking, eg. 'left_x_data'
%   EyePosYName             - the name of the channel containing the
%                             y-coordinate of the eye-tracking, eg. 'left_y_data'
%   trailColor              - colormap or 2 x 3 matrix of colors.
%                             Takes also standardized maps, eg. 'summer', 'magma', or
%                             'cividis' (DEFAULT = 'magma') (OPTIONAL ARGUMENT)
%   names                   - provide the names as cells for the motions to appear
%                             in the plot (eg. {'Incorrect'; 'Correct'} (OPTIONAL ARGUMENT)
%   norm                    - normalize the data so that it fits between [-1 1] (DEFAULT = 1)
%   figpos                  - set the figure size and position (OPTIONAL ARGUMENT)
%
% Outputs:
%   [coorddata, handleList] - eye tracking data structure with fieldnames,
%                             and the figure handles when the figures are open.
%
% Example:
%   [coorddata, handleList]= bemobil_plot_eyetracking({ALLEEG(1) ALLEEG(2)},...
%                               trailColor, 'magma', ...
%                               'EyePosXName', 'CombinedEyes_X',
%                               'EyePosYName', 'CombinedEyes_Y',
%                               norm, 1)
%
% See also:
%   EEGLAB, bemobil_plot_trail, bemobil_plot_heatmap, bemobil_plot_motion
%
% Author: Zak Djebbara, April 2022

function [coorddata, handleList] = bemobil_plot_eyetracking(EEG, varargin)

if nargin == 0; help bemobil_plot_eyetracking; return; end

p = inputParser;
epoched = [];
handleList = [];

addRequired(p, 'EEG');
addParameter(p, 'names','');
addParameter(p, 'figPos', [200 200 800 300], @isnumeric);
addParameter(p, 'trailColor', 'magma');
addParameter(p, 'EyePosXName', 'Eyetracking_Screen_X_both', @ischar);
addParameter(p, 'EyePosYName', 'Eyetracking_Screen_Y_both', @ischar);
addParameter(p, 'norm', 1, @isnumeric);

parse(p, EEG, varargin{:});
figPos      = p.Results.figPos;
trailColor = p.Results.trailColor;
EyePosXName = p.Results.EyePosXName;
EyePosYName = p.Results.EyePosYName;
names       = p.Results.names;
norm        = p.Results.norm;

% Creating names if not provided
if isempty(names)
    if iscell(EEG)
        for i1 = 1:size(EEG,2)
            names{i1} = ['Channel ' num2str(i1)];
        end
    elseif ~iscell(names)
        error('Please make sure you input the EEG struct as a cell: {EEG} or {ALLEEG}');
    end
elseif ~isempty(names)
    if iscell(names)
        if ~ischar(names{1})
            error('Please use characters to for the names');
        end
    elseif ~iscell(names)
        error(['Please place the names in cells, eg.:' newline ...
            'names{1} = "somename"; names{2} = "someothername";']);
    end
end

% Setting the colormap
% Setting up the colormap
    try
        if ischar(trailColor)
            cmap = eval([trailColor '(2);']);
            xcmap = cmap(1,:);
            ycmap = cmap(2,:);
        elseif isnumeric(trailColor)
            cmap = trailColor;
            xcmap = cmap(1,:);
            ycmap = cmap(2,:);
        end
    catch
        error(['Colormap not recognized. Please use standard maps, eg. "magma", "summer", or "winter"'...
            newline 'or matrix of at least two rows, eg. [0 0 0; .5 .5 .5];']);
    end
figure;
t = tiledlayout(2, length(EEG));
t.TileSpacing = 'compact';
t.Padding = 'compact';
for i = 1:length(EEG)
    % Checking if data is epoched or not
    if isempty(EEG{i}.epoch)
        warning(['No epochs detected. Generates continuous plot.']);
        epoched = 0;
    elseif ~isempty(EEG{i}.epoch)
        warning(['Epoched data detected. Generates average trail of data.']);
        epoched = 1;
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
        if norm; xdata = normalize(xdata, 'range', [-1 1]); xdata = normalize(xdata, 'center');end
        xdatadfdt   = normalize(gradient(xdata),'range',[min(xdata) max(xdata)]);
        ydata       = mean(EEG{i}.data,3);
        ydata       = ydata(ychanidx,:);
        if norm; ydata = normalize(ydata, 'range', [-1 1]); ydata = normalize(ydata, 'center');end
        ydatadfdt   = normalize(gradient(ydata),'range',[min(ydata) max(ydata)]);
        tickstep    = length(EEG{i}.times)/5;
        ticksX      = 1:round(tickstep):length(xdata);
        ticksX      = [ticksX, length(xdata)];
    elseif ~epoched
        xdata       = EEG{i}.data(xchanidx,:);
        if norm; xdata = normalize(xdata, 'range', [-1 1]); xdata = normalize(xdata, 'center');end
        xdatadfdt   = normalize(gradient(xdata),'range',[min(xdata) max(xdata)]);
        ydata       = EEG{i}.data(ychanidx,:);
        if norm; ydata = normalize(ydata, 'range', [-1 1]); ydata = normalize(ydata, 'center');end
        ydatadfdt   = normalize(gradient(ydata),'range',[min(ydata) max(ydata)]);
        
        tickstep    = length(EEG{i}.times)/5;
        ticksX      = 1:round(tickstep):length(xdata);
        ticksX      = [ticksX, length(xdata)];
    end

    % Plotting
    h = nexttile;
        if norm; xdata = normalize(xdata, 'range', [-1 1]);end
        plot(xdata,'Color',xcmap);
        hold on;
        plot(xdatadfdt,':','Color',xcmap);
        if epoched; xline([find(EEG{i}.times == 0)],':','Onset');end
        set(h,'XTickLabel',(num2cell(EEG{i}.times(ticksX))));
        set(gca, 'XTick',ticksX);
        title([names{i} newline 'X-direction.']);
        xlabel('time (ms)'); ylabel('x coordinate');
        xlim([0 length(xdata)]);
        handleList = [handleList; gcf];
        hold off;
        
    g = nexttile;
        plot(ydata,'Color',ycmap);
        hold on;
        plot(ydatadfdt,':','Color',ycmap);
        if epoched; xline([find(EEG{i}.times == 0)],':','Onset');end
        set(g,'XTickLabel',(num2cell(EEG{i}.times(ticksX))));
        set(gca, 'XTick',ticksX);
        title([names{i} newline 'Y-direction.']);
        xlabel('time (ms)'); ylabel('y coordinate');
        xlim([0 length(ydata)]);
        handleList = [handleList; gcf];
        hold off;
        
        % Saving the data
        coorddata{i}.xdata=xdata;
        coorddata{i}.ydata=ydata;
end
set(gcf, 'Color', 'w', 'Position', figPos);
end