% bemobil_plot_heatmap()- Plots the trail map of the selected eye tracking
%                         data in the EEG data strucutre. Takes both
%                         continuous and epoched data. The continuous data
%                         is not recommended as it may take a long time to
%                         display. The thickness of the trail map
%                         represents the velocity of the movement (thin =
%                         slow, thick = fast).
%
% Usage:
%   >> [trailmap, handleList] = bemobil_plot_heatmap(EEG, trailFactor, trailColor, trailSize,...
%                               EyePosXName, EyePosYName, norm, figPos)
% Inputs:
%   EEG                     - current EEGLAB EEG structure. Takes only a
%                             single EEG struct at the time.
%   colorMap                - colormap or N x 3 matrix of colors, where N > 2.
%                             Takes also standardized maps, eg. 'summer', 'magma', or
%                             'cividis'. Can also be set to each input EEG dataset
%                             in the following way; {'magma' 'cividis'}. This sets
%                             the first as 'magma, the other as 'cividis'.
%                             (DEFAULT = 'magma') (OPTIONAL ARGUMENT)
%   EyePosXName             - the name of the channel containing the
%                             x-coordinate of the eye-tracking, eg. 'left_x_data'
%   EyePosYName             - the name of the channel containing the
%                             y-coordinate of the eye-tracking, eg. 'left_y_data'
%   addAxes                 - Adds cross-axes over the heatmap to better
%                             visualize what's up/down and right/left (OPTIONAL ARGUMENT)
%   norm                    - normalize the data so that it fits between [-1 1] (DEFAULT = 1)
%   figpos                  - set the figure size and position (OPTIONAL ARGUMENT)
%
% Outputs:
%   [heatmapdata, handleList]  - heatmap data structure with fieldnames,
%                             and the figure handles when the figures are open.
%
% Example:
%   [heatmapdata, handleList]   = bemobil_plot_heatmap({ALLEEG(1) ALLEEG(2)},
%                               'EyePosXName', 'CombinedEyes_X',
%                               'EyePosYName', 'CombinedEyes_Y',
%                               'colorMap', {'viridis' 'magma'},
%                               'trailsize', [1 50], 'norm', 1);
%
% See also:
%   EEGLAB, bemobil_plot_eyetracking, bemobil_plot_trail, bemobil_plot_motion
%
% Author: Zak Djebbara, April 2022

function [heatmapdata, handleList] = bemobil_plot_heatmap(EEG,varargin)

if nargin == 0; help bemobil_heatmap; return; end

p = inputParser;
epoched = [];
handleList = [];

addRequired(p, 'EEG');

addParameter(p, 'EyePosXName', 'Eyetracking_Screen_X_both', @ischar);
addParameter(p, 'EyePosYName', 'Eyetracking_Screen_Y_both', @ischar);
addParameter(p, 'names', '');
addParameter(p, 'colorMap', 'magma');
addParameter(p, 'addAxes', 1), @isnumeric;
addParameter(p, 'figPos', [300 300 600 400]), @isnumeric;
addParameter(p, 'norm', 1, @isnumeric);

parse(p, EEG, varargin{:});
EyePosXName = p.Results.EyePosXName;
EyePosYName = p.Results.EyePosYName;
names       = p.Results.names;
colorMap    = p.Results.colorMap;
addAxes     = p.Results.addAxes;
figPos      = p.Results.figPos;
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

% Checking if data is epoched or not
if isempty(EEG{1}.epoch)
    warning(['No epochs detected. Generates heatmap for continuous data.'...
        newline 'This is not adviced. Please be patient as this may take some time. Ctrl+C to cancel operation.']);
    epoched = 0; plottitle = "continuous data";
elseif ~isempty(EEG{1}.epoch)
    warning(['Epoched data detected. Generates average heatmap of data.']);
    epoched = 1; plottitle = "epoched data";
end
pause(.1)

figure('color','w');
tiledlayout('flow', 'TileSpacing','compact' , 'Padding', 'compact');
for i = 1:length(EEG)
    
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
        data = mean(EEG{i}.data,3);
        xdata = data(xchanidx,:);
        ydata = data(ychanidx,:);
    elseif ~epoched
        xdata = EEG{i}.data(xchanidx,:);
        ydata = EEG{i}.data(ychanidx,:);
    end

    % Plotting 
    gridx1 = -1:.01:1;
    gridx2 = -1:.01:1;
    [x1,x2] = meshgrid(gridx1, gridx2);
    x1 = x1(:);
    x2 = x2(:);
    xi = [x1 x2];
    if norm; data = [normalize(xdata,'range',[-1 1]); normalize(ydata,'range',[-1 1])];end
    if ~norm; data = [xdata; ydata];end
    
    % Generating heatmap data
    nexttile(i);
    [f,x]=ksdensity(data', xi);
    heatmapdata{i}.probabilitydensityestimate = reshape(f,size(meshgrid(gridx1, gridx2)));
    heatmapdata{i}.xi = x;
    
    % Setting the colormap
    if iscell(colorMap)
        if size(colorMap,2) < size(EEG,2)
            colorMap = repmat(colorMap, size(EEG));
        elseif size(colorMap,2) > size(EEG,2)
            error('There are more colors than EEG data provided. Please make sure they match.');
        end
    end
    
    try
        if iscell(colorMap)
            if ischar(colorMap{1})
                cmap = eval([colorMap{i} '(' num2str(length(x(:,1))) ');']);
            elseif isnumeric(colorMap{i})
                cmap = bemobil_makecmap(colorMap{i}, length(x(:,1)));
            end
        elseif ~iscell(colorMap)
            if ischar(colorMap)
                cmap = eval([colorMap '(' num2str(length(x(:,1))) ');']);
            elseif isnumeric(colorMap)
                cmap = bemobil_makecmap(colorMap, length(x(:,1)));
            end
        end
    catch
        error(['Colormap not recognized. Please use standardized colormaps or' ... 
            newline 'a N x 3, where N > 2, matrix']);
    end
    
    imagesc(x(:,1),x(:,2),reshape(f,size(meshgrid(gridx1, gridx2))));
    set(gca,'YDir','normal','ColorMap',cmap)
    set(gcf, 'Position', figPos);
    title(['Heatmap (KDE) of ' names{i} '.']);
    shading interp; xlabel('x coordinate'); ylabel('y coordinate');
    view(0,90);
    axis('image');
    if addAxes
        hold on
        plot(ylim,[0 0],'w');
        plot([0 0],xlim,'w');
    end
    handleList = [handleList; gcf];
end
set(gcf, 'Color', 'w');
end