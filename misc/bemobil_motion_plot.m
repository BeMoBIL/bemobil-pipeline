% Some help text
% motionChannels = {['motion_channel_X']; ['motion_channel_Y']; ['motion_channel_Z']};
%
% Zak, 2022
function [motiondata] = bemobil_motion_plot(EEG, motionChannels, varargin)
if nargin == 0; help bemobil_motion_plot; end

p = inputParser;
epoched = [];

addRequired(p, 'EEG');
addRequired(p, 'motionChannels', @iscellstr);

addParameter(p, 'TwoDcolorMap', 'magma');
addParameter(p, 'ThreeDcolorMap', 'magma');
addParameter(p, 'figPos', [200 200 1000 400]), @isnumeric;
addParameter(p, 'axeson', 0, @isnumeric);

parse(p, EEG, motionChannels, varargin{:});
TwoDcolorMap    = p.Results.TwoDcolorMap;
ThreeDcolorMap  = p.Results.ThreeDcolorMap;
figPos          = p.Results.figPos;
axeson          = p.Results.axeson;

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
for i1 = 1:length(EEG.chanlocs)
    if strcmp(EEG.chanlocs(i1).labels, motionChannels{1}); xchanidx = i1;end
    if strcmp(EEG.chanlocs(i1).labels, motionChannels{2}); ychanidx = i1;end
    if strcmp(EEG.chanlocs(i1).labels, motionChannels{3}); zchanidx = i1;end
end
if size(ychanidx,2)> 1 || size(xchanidx,2)> 1 || size(zchanidx,2) > 1
    error('More than one channel share the same name for motion channel(s).');
end
if isempty(ychanidx) || isempty(xchanidx); error('Motion channel(s) could not be found.');end

% Identifying data
if epoched
    data = mean(EEG.data,3);
    xdata = data(xchanidx,:);
    ydata = data(ychanidx,:);
    zdata = data(zchanidx,:);

    xdatadfdt = normalize(gradient(xdata),'range',[min(xdata) max(xdata)]);
    ydatadfdt = normalize(gradient(ydata),'range',[min(ydata) max(ydata)]);
    zdatadfdt = normalize(gradient(zdata),'range',[min(zdata) max(zdata)]);
elseif ~epoched
    xdata = EEG.data(xchanidx,:);
    ydata = EEG.data(ychanidx,:);
    zdata = EEG.data(zchanidx,:);
    
    xdatadfdt = normalize(gradient(xdata),'range',[min(xdata) max(xdata)]);
    ydatadfdt = normalize(gradient(ydata),'range',[min(ydata) max(ydata)]);
    zdatadfdt = normalize(gradient(zdata),'range',[min(zdata) max(zdata)]);
end

% Plotting
% Setting the colormap
try
    cmap2 = eval([TwoDcolorMap '(' num2str(length(motionChannels)*2) ');']);
catch
    if size(TwoDcolorMap) == [size(motionChannels,1)*2 3]
        cmap2 = TwoDcolorMap;
    else
        error(['Colormap not recognized. Please use [number of motionchannels*2 x 3].'...
            newline 'Or simply "magma" or "cividis".']);
    end
end

% Determine the xticks
tickstep = length(EEG.times)/5;
ticksX = 1:round(tickstep):length(xdata);
ticksX = [ticksX, length(xdata)];
if epoched; ticksXlabels = (num2cell(EEG.times(ticksX))); labelsX = 'ms'; end
if ~epoched; ticksXlabels = (num2cell(EEG.times(ticksX)/1000)); labelsX = 's'; end


% Design plot
t = tiledlayout(3,2);

% X coordinate
nexttile(1);
plot(xdata, 'Color', cmap2(1,:));
hold on;
plot(xdatadfdt, ':', 'Color', cmap2(2,:));
if epoched; xline(find(EEG.times == 0),':','Onset');end
xticklabels(ticksXlabels);
xticks(ticksX);
xlabel(['time (' labelsX ')']); 
if isfield(EEG.chanlocs(xchanidx), 'units') || isfield(EEG.chanlocs(xchanidx), 'unit')
    yunit = EEG.chanlocs(xchanidx).units;
end
ylabel(['distance (' yunit ')']);
title(['X-direction of ' plottitle ' motion data.']);

% Y coordinate
nexttile(3);
plot(ydata, 'Color', cmap2(3,:));
hold on;
plot(ydatadfdt, ':', 'Color', cmap2(4,:));
if epoched; xline(find(EEG.times == 0),':','Onset');end
xticklabels(ticksXlabels);
xticks(ticksX);
xlabel(['time (' labelsX ')']); 
if isfield(EEG.chanlocs(xchanidx), 'units') || isfield(EEG.chanlocs(xchanidx), 'unit')
    yunit = EEG.chanlocs(xchanidx).units;
end
ylabel(['distance (' yunit ')']);
title(['Y-direction of ' plottitle ' motion data.']);

% Z coordinate
nexttile(5);
plot(zdata, 'Color', cmap2(5,:));
hold on;
plot(zdatadfdt, ':', 'Color', cmap2(6,:));
if epoched; xline(find(EEG.times == 0),':','Onset');end
xticklabels(ticksXlabels);
xticks(ticksX);
xlabel(['time (' labelsX ')']); 
if isfield(EEG.chanlocs(xchanidx), 'units') || isfield(EEG.chanlocs(xchanidx), 'unit')
    yunit = EEG.chanlocs(xchanidx).units;
end
ylabel(['distance (' yunit ')']);
title(['Z-direction of ' plottitle ' motion data.']);

% 3D trajectory plot
% Setting the colormap
try
    cmap3 = eval([ThreeDcolorMap '(' num2str(length(xdata)) ');']);
catch
    if size(ThreeDcolorMap) == [length(xdata) 3]
        cmap3 = ThreeDcolorMap;
    else
        error(['Colormap not recognized. Please use [length(data) x 3].'...
            newline 'Or simply "summer", "winter", "magma" or "cividis".']);
    end
end

% Determining the size of the trajectory
ThreeDdata = [xdata;ydata;zdata];
datadfdt = gradient(ThreeDdata);
avgdatadfdt = mean(datadfdt);
sz = normalize(avgdatadfdt, 'range', [1 50]);

% Plotting
nexttile(2, [3 1]);
scatter3(xdata,ydata,zdata, sz, cmap3, 'filled'); colormap(cmap3);
cb = colorbar('XTickLabel',{'Start','Midway','End'},'XTick',[0 .5 1], 'Location', 'eastoutside');
cb.Label.String = 'temporal development';
% set(cb, 'YDir', 'reverse');
title(['Trajectory of ' plottitle ' motion data.']);

% Finalize plot
title(t, ['BeMoBIL motion data plot of:' newline...
    EEG.chanlocs(xchanidx).labels ', ' EEG.chanlocs(xchanidx).labels ', and ' EEG.chanlocs(zchanidx).labels],...
    'Interpreter', 'none', 'fontweight', 'bold', 'fontsize',9);
set(gcf, 'Position', figPos);
set(gca, 'Color', 'none');
axis('image');
if axeson; axis on;end
if ~axeson; axis off;end

% Outputs
motiondata.xdata = xdata;
motiondata.ydata = ydata;
motiondata.zdata = zdata;

motiondata.xdatadfdt = xdatadfdt;
motiondata.ydatadfdt = ydatadfdt;
motiondata.zdatadfdt = zdatadfdt;

motiondata.averagedfdt = avgdatadfdt;
end