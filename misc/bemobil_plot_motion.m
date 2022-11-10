% bemobil_plot_motion() - Plots the selected motion data in the EEG
%                         data structure. Takes both continuous and epoched
%                         datasets. The thickness of the 3D plot represents
%                         the velocity of the movement (thin = slow, thick
%                         = fast).
%
% Usage:
%   >> [motiondata, fHandle]=   bemobil_plot_motion( EEG , motionChannels , names , range ,...
%                               TwoDcolorMap , ThreeDcolorMap , axeson , figPos)
% Inputs:
%   EEG                     - current EEGLAB EEG structure. Takes only a
%                             single EEG struct at the time.
%   motionChannels          - name of the X, Y, and Z motion channels as provided in
%                             the EEG.chanlocs.labels'-structure (e.g. {'motion_channel_X';
%                             'motion_channel_Y'; 'motion_channel_Z'})
%   names                   - provide the names as cells for the motions to appear
%                             in the plot (eg. {'right hand'} or {'right hand'; 'head'}
%   range                   - set the range for the thickness of 3D (e.g. [1 30], where [min max])
%                             trajectory plot. The thickness represents the velocity of the movement.
%   TwoDcolorMap            - (2 x 3) colormap or standard colormaps (e.g. 'summer')
%                             for the timeseries plots (e.g. [1 1 1; 0 0 0]) (DEFAULT = 'magma')
%   ThreeDcolorMap          - (2 x 3) colormap or standard colormaps (e.g. 'summer')
%                             for the 3D trajectory plot (e.g. [1 1 1; 0 0 0]) (DEFAULT = 'magma')
%   axeson                  - set whether to plot axes on 3D plot (OPTIONAL ARGUMENT)
%   figpos                  - set the figure size and position (OPTIONAL ARGUMENT)
%
% Outputs:
%   motiondata              - motiondata structure with fieldnames
%
% Example:
%   motiondata = bemobil_plot_motion(EEG , {'right_hand_x'; 'right_hand_y'; 'right_hand_z'},...
%                axeson , 0 , 'TwoDcolorMap' , 'magma' , 'ThreeDcolorMap' , 'viridis')
%
% See also:
%   EEGLAB, bemobil_plot_eyetracking, bemobil_plot_heatmap, bemobil_plot_patterns
%
% Author: Zak Djebbara, April 2022
function [motiondata, fHandle] = bemobil_plot_motion(EEG, motionChannels, varargin)
if nargin == 0; help bemobil_plot_motion; return; end

p = inputParser;
epoched = [];

addRequired(p, 'EEG');
addRequired(p, 'motionChannels', @iscell);

addParameter(p, 'names', []);
addParameter(p, 'range', [1 20], @isnumeric);
addParameter(p, 'TwoDcolorMap', 'magma');
addParameter(p, 'ThreeDcolorMap', 'magma');
addParameter(p, 'figPos', [200 200 1100 400]), @isnumeric;
addParameter(p, 'axeson', 1, @isnumeric);

parse(p, EEG, motionChannels, varargin{:});
names           = p.Results.names;
range           = p.Results.range;
TwoDcolorMap    = p.Results.TwoDcolorMap;
ThreeDcolorMap  = p.Results.ThreeDcolorMap;
figPos          = p.Results.figPos;
axeson          = p.Results.axeson;


% Making sure 3 coordinates are provided
errormsg = 'Provide XYZ coordinates for each motion channel.';
if iscell(motionChannels)
    for i1 = 1:length(motionChannels)
        if length(motionChannels{i1}) < 3; error(errormsg); end
        if ~ischar(motionChannels{i1}{1}); error('Please use characters for the motion channels.'); end
    end
else
    if size(motionChannels,1) < 3; error(errormsg); end
end

% Creating names if not provided
if isempty(names)
    if iscell(motionChannels)
        for i1 = 1:size(motionChannels,2)
            names{i1} = ['Channel ' num2str(i1)];
        end
    elseif ~iscell(motionChannels)
        names{1} = [names; 'Channel 1'];
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

% Checking if data is epoched or not
if isempty(EEG.epoch)
    warning(['No epochs detected. Generates motion plot for continuous data.']);
    epoched = 0; plottitle = 'continuous';
elseif ~isempty(EEG.epoch)
    warning(['Epoched data detected. Generates averaged motion plot of data.']);
    epoched = 1; plottitle = 'averaged and epoched';
end
pause(.1)

% Identifying the channel names
xchanidx = [];
ychanidx = [];
zchanidx = [];
for i1 = 1:length(motionChannels)
    for i2 = 1:length(EEG.chanlocs)
        if strcmp(EEG.chanlocs(i2).labels, motionChannels{i1}{1}); xchanidx{i1}{1} = i2;end
        if strcmp(EEG.chanlocs(i2).labels, motionChannels{i1}{2}); ychanidx{i1}{1} = i2;end
        if strcmp(EEG.chanlocs(i2).labels, motionChannels{i1}{3}); zchanidx{i1}{1} = i2;end
    end
    if size(ychanidx{i1},2)> 1 || size(xchanidx{i1},2)> 1 || size(zchanidx{i1},2) > 1
        error('More than one channel share the same name for motion channel(s).');
    end
    if isempty(ychanidx) || isempty(xchanidx) || isempty(zchanidx); error('Motion channel(s) could not be found.');end
end

% Identifying data
for i1 = 1:length(motionChannels)
    if epoched
        data = mean(EEG.data,3);
        xdata{i1} = data(xchanidx{i1}{1},:);
        ydata{i1} = data(ychanidx{i1}{1},:);
        zdata{i1} = data(zchanidx{i1}{1},:);

        xdatadfdt{i1} = normalize(gradient(xdata{i1}),'range',[min(xdata{i1}) max(xdata{i1})]);
        ydatadfdt{i1} = normalize(gradient(ydata{i1}),'range',[min(ydata{i1}) max(ydata{i1})]);
        zdatadfdt{i1} = normalize(gradient(zdata{i1}),'range',[min(zdata{i1}) max(zdata{i1})]);
    elseif ~epoched
        xdata{i1} = EEG.data(xchanidx{i1}{1},:);
        ydata{i1} = EEG.data(ychanidx{i1}{1},:);
        zdata{i1} = EEG.data(zchanidx{i1}{1},:);

        xdatadfdt{i1} = normalize(gradient(xdata{i1}),'range',[min(xdata{i1}) max(xdata{i1})]);
        ydatadfdt{i1} = normalize(gradient(ydata{i1}),'range',[min(ydata{i1}) max(ydata{i1})]);
        zdatadfdt{i1} = normalize(gradient(zdata{i1}),'range',[min(zdata{i1}) max(zdata{i1})]);
    end
end

% Plotting
% Setting the colormap
try
    if iscell(TwoDcolorMap)
        if ischar(TwoDcolorMap{1})
            for i2 = 1:size(motionChannels,2)
                cmap2{i2} = eval([TwoDcolorMap{1} '(' num2str(length(motionChannels{i2})*2) ');']);
            end
        else 
            try
                cmap2{1} = bemobil_makecmap(TwoDcolorMap{:}, 6);
            catch
                error('Cannot recognize the colormap. If in doubt, please simply use "summer", "winter," or "magma."');
            end
        end
    elseif ~iscell(TwoDcolorMap)
        cmap2{1} = eval([TwoDcolorMap '(6);']);
    end
    
catch
    % in case of multiple inputs
    if iscell(TwoDcolorMap) && iscell(motionChannels) 
        for i2 = 1:length(motionChannels)
            try
                if size(TwoDcolorMap{i2},1) == size(motionChannels{i2},1)*2
                    if i2 == length(TwoDcolorMap); cmap2 = TwoDcolorMap; end
                else
                    error(['The provided colormap does not match the number of motions channels provided.'...
                        newline 'If in doubt, please simply use "summer", "winter," or "magma."']);
                end
            catch
                error(['The provided colormap does not match the number of motions channels provided.'...
                    newline 'If in doubt, please simply use "summer", "winter," or "magma."']);
            end
        end
        
    % in case of singular input
    elseif ~iscell(TwoDcolorMap)
        if ~ischar(TwoDcolorMap)
            try
                cmap2{1} = bemobil_makecmap(TwoDcolorMap, 6);
            catch
                error(['Non-standard colormap detected. Please use at least two colors as input, e.g.'...
                    newline '[0 0 0; .5 .5 .5]; for an automatic generation of colors.']);
            end
        end
    elseif ~iscell(TwoDcolorMap) && iscell(motionChannels) || iscell(TwoDcolorMap) && ~iscell(motionChannels)
        error(['The provided colormap and motion channels do not match in structure.'...
            newline 'Make sure they are either both cells or arrays']);
    end
end

% Determine the xticks
tickstep = length(EEG.times)/5;
ticksX = 1:round(tickstep):length(xdata{1});
ticksX = [ticksX, length(xdata{1})];
if epoched; ticksXlabels = (num2cell(EEG.times(ticksX))); labelsX = 'ms'; end
if ~epoched; ticksXlabels = (num2cell(EEG.times(ticksX)/1000)); labelsX = 's'; end


% Design plot
figure;
fHandle = tiledlayout(3,2, 'TileSpacing','compact' , 'Padding', 'compact');
if iscell(motionChannels); motionsize = size(motionChannels,2); end
if ~iscell(motionChannels); motionsize = 1; end

% X coordinate
nexttile(1);
plots = [];
for i3 = 1:motionsize
    a = plot(xdata{i3}, 'Color', cmap2{i3}(1,:), 'DisplayName', names{i3}, 'LineWidth', 1.5);
    hold on;
    b = plot(xdatadfdt{i3}, ':', 'Color', cmap2{i3}(2,:), 'DisplayName', [names{i3} 'df/dt'], 'LineWidth', 1.5);
    if epoched; xline(find(EEG.times == 0),'-','LineWidth', 1);end
    xticklabels(ticksXlabels);
    xticks(ticksX);
    xlabel(['time (' labelsX ')']);
    if isfield(EEG.chanlocs(xchanidx{i3}{1}), 'units') || isfield(EEG.chanlocs(xchanidx{i3}{1}), 'unit')
        yunit = EEG.chanlocs(xchanidx{i3}{1}).units;
    else
        yunit = 'unknown';
    end
    ylabel(['distance (' yunit ')']);
    title(['X-direction of ' plottitle ' motion data.']);
    plots = [plots; a;b];
    legend(plots, 'location', 'westoutside');
end

% Y coordinate
nexttile(3);
plots = [];
for i3 = 1:motionsize
    a = plot(ydata{i3}, 'Color', cmap2{i3}(3,:), 'DisplayName',names{i3}, 'LineWidth', 1.5);
    hold on;
    b = plot(ydatadfdt{i3}, ':', 'Color', cmap2{i3}(4,:), 'DisplayName', [names{i3} 'df/dt'], 'LineWidth', 1.5);
    if epoched; xline(find(EEG.times == 0),'-','LineWidth', 1);end
    xticklabels(ticksXlabels);
    xticks(ticksX);
    xlabel(['time (' labelsX ')']); 
    if isfield(EEG.chanlocs(xchanidx{i3}{1}), 'units') || isfield(EEG.chanlocs(xchanidx{i3}{1}), 'unit')
        yunit = EEG.chanlocs(xchanidx{i3}{1}).units;
        else
        yunit = 'unknown';
    end
    ylabel(['distance (' yunit ')']);
    title(['Y-direction of ' plottitle ' motion data.']);
    plots = [plots; a;b];
    legend(plots, 'location', 'westoutside');
end

% Z coordinate
nexttile(5);
plots = [];
for i3 = 1:motionsize
    a = plot(zdata{i3}, 'Color', cmap2{i3}(5,:), 'DisplayName',names{i3}, 'LineWidth', 1.5);
    hold on;
    b = plot(zdatadfdt{i3}, ':', 'Color', cmap2{i3}(6,:), 'DisplayName', [names{i3} 'df/dt'], 'LineWidth', 1.5);
    if epoched; xline(find(EEG.times == 0),'-','LineWidth', 1);end
    xticklabels(ticksXlabels);
    xticks(ticksX);
    xlabel(['time (' labelsX ')']); 
    if isfield(EEG.chanlocs(xchanidx{i3}{1}), 'units') || isfield(EEG.chanlocs(xchanidx{i3}{1}), 'unit')
        yunit = EEG.chanlocs(xchanidx{i3}{1}).units;
    else
        yunit = 'unknown';
    end
    ylabel(['distance (' yunit ')']);
    title(['Z-direction of ' plottitle ' motion data.']);
    plots = [plots; a;b];
    legend(plots, 'location', 'westoutside');
end

% 3D trajectory plot
% Setting the colormap
try
    if iscell(ThreeDcolorMap)
        try
            for i2 = 1:size(motionChannels,2)
                cmap3{i2} = eval([ThreeDcolorMap{i2} '(' num2str(length(xdata{i2})) ');']);
            end
        catch
            cmap3{1} = bemobil_makecmap(ThreeDcolorMap{:}, length(xdata{1}));
        end
    elseif ~iscell(ThreeDcolorMap)
        if ischar(ThreeDcolorMap)
            cmap3{1} = eval([ThreeDcolorMap{1} '(' num2str(length(xdata{1})) ');']);
        else
            try
                cmap3{1} = bemobile_makecmap(ThreeDcolorMap, length(xdata{1}));
            catch
                error(['Non-standard colormap detected. Please use at least two colors as input, e.g.'...
                    newline '[0 0 0; .5 .5 .5]; for an automatic generation of colors.']);
            end
        end
    end
    
catch
    if ischar(ThreeDcolorMap)
        cmap3{1} = eval([ThreeDcolorMap '(' num2str(length(xdata{1})) ');']);
    else
        try
            cmap3{1} = bemobil_makecmap(ThreeDcolorMap, length(xdata{1}));
        catch 
            error(['Colormap not recognized. Please use [length(data) x 3].'...
                newline 'Or simply "summer", "winter", "magma" or "cividis".']);
        end
    end
end

% Determining the size of the trajectory
for i3 = 1:motionsize
    ThreeDdata{i3} = [xdata{i3};ydata{i3};zdata{i3}];
    datadfdt{i3} = gradient(ThreeDdata{i3});
    avgdatadfdt{i3} = mean(datadfdt{i3});
    sz{i3} = normalize(avgdatadfdt{i3}, 'range', range);
end

% Plotting
nexttile(2, [3 1]);
plots = [];
axeslist = [];
for i3 = 1:motionsize
    a = scatter3(xdata{i3},ydata{i3},zdata{i3}, sz{i3}, cmap3{i3}, 'filled'); colormap(cmap3{i3});
    hold on;
    title(['Trajectory of ' names{i3} ' motion data.']);
    plots = [plots; a;];
%     axeslist = [axeslist; ax1;];
%     cb = colorbar(gca, 'Location', 'southoutside');
%     legend(plots, 'location', 'eastoutside');
    cb = colorbar(gca, 'Location', 'southoutside','XTickLabel',{'Start','Midway','End'},'XTick',[0 .5 1]);
%     cb.Label.String = ['time of ' names{i3}];
end

% Finalize plot
title(fHandle, 'BeMoBIL motion data plot', 'Interpreter', 'none', 'fontweight', 'bold', 'fontsize',9);
set(gcf, 'Position', figPos);
set(gca, 'Color', 'none');
axis('image');
if axeson; axis on;end
if ~axeson; axis off;end

% Outputs
for i3 = 1:motionsize
    motiondata.xdata = xdata{i3};
    motiondata.ydata = ydata{i3};
    motiondata.zdata = zdata{i3};

    motiondata.xdatadfdt = xdatadfdt{i3};
    motiondata.ydatadfdt = ydatadfdt{i3};
    motiondata.zdatadfdt = zdatadfdt{i3};

    motiondata.averagedfdt = avgdatadfdt{i3};
end
set(gcf, 'Color', 'w');
end