% bemobil_detect_motion_startstops detects motion starts and stops based on a coarse and fine threshold of one or more
% given channels. If more than one channel is given, the square root of the sum of squares is taken as the detection
% data. A movement is detected firstly based on a coarse threshold of a given quantile of the data, then a fine
% threshold is applied afterwards based on the fine threshold quantile in a buffer around the detected coarse movement.
% Events and used parameters are stored in the EEG.etc struct!
%
% The algorithm was first used and published in Gramann, K., Hohlefeld, F.U., Gehrke, L. et al. Human cortical dynamics
% during full-body heading changes. Sci Rep 11, 18186 (2021). https://doi.org/10.1038/s41598-021-97749-8
%
% Input:
%       EEG_motion                  - EEGLAB dataset containing motion data
%       idx_detect                  - OPTIONAL channel index to detect (default = 1)
%       eventlabel                  - OPTIONAL label the events should receive (plus the suffixes "start" and "stop") (default = motion)
%       movement_threshold          - OPTIONAL coarse movement threshold. This is the quantile of the channel that needs 
%                                       to be exceeded for a motion to be detected (default = 0.65)
%       movement_threshold_fine     - OPTIONAL fine movement threshold that is used when the coarse detection is positive (default = 0.05)
%       min_duration                - OPTIONAL minimal movement duration in seconds (default = 0)
%       detection_buffer            - OPTIONAL buffer that is used when a coarse movement is detected (default = 2)
%       search_timerange            - OPTIONAL timerange in which the search should happen (useful if different
%                                       conditions exist, default = all data)
%       createplots                 - OPTIONAL boolean whether plots should be created (default = 1)
%
% Output:
%       EEG_motion          - EEGLAB dataset containing detected events
%       plothandles         - handle to the created plots to allow script-based saving and closing
%
% Example: 
%         [EEG_motion, plothandles] = bemobil_detect_motion_startstops(EEG_motion,idx_detect,eventlabel,movement_threshold,...
%             movement_threshold_fine,min_duration,detection_buffer,search_timerange,createplots)
% 
% Authors: Marius Klug, 2022

function [EEG_motion, plothandles] = bemobil_detect_motion_startstops(EEG_motion,idx_detect,eventlabel,movement_threshold,...
    movement_threshold_fine,min_duration,detection_buffer,search_timerange,createplots)

if nargin == 0
    help bemobil_ECG_analysis
    return
end

if ~isfield(EEG_motion,'data')
    error('Wrong data type for EEG_motion input! Needs to be an EEGLAB struct.')
end

if ~exist('idx_detect','var') || ~isvector(idx_detect)
    disp('Assuming only one channel exists in the data.')
    idx_detect = 1;
end

if ~exist('eventlabel','var') || ~ischar(eventlabel)
    eventlabel = 'motion';
    disp(['Using default event label "' eventlabel '"'])
end

if ~exist('movement_threshold','var') || ~isscalar(movement_threshold)
    movement_threshold = 0.65;
    disp(['Using the default movement detector threshold of ' num2str(movement_threshold) '.'])
end

if ~exist('movement_threshold_fine','var') || ~isscalar(movement_threshold_fine)
    movement_threshold_fine = 0.05;
    disp(['Using the default fine movement detector threshold of ' num2str(movement_threshold_fine) '.'])
end

if ~exist('min_duration','var') || ~isscalar(min_duration)
    min_duration = 0;
    disp(['Using the default minimum motion duration of ' num2str(min_duration) '.'])
end

if ~exist('detection_buffer','var') || ~isscalar(detection_buffer)
    detection_buffer = 2;
    disp(['Using the default fine detection buffer of ' num2str(detection_buffer) '.'])
end

if ~exist('search_timerange','var') || ~isvector(search_timerange)
    search_timerange = [1 EEG_motion.pnts];
    disp('Using the default search time range of the entire data set.')
end

if ~exist('createplots','var') || ~isscalar(createplots)
    createplots = 1;
end

detection_buffer_use = detection_buffer*EEG_motion.srate;

%% apply motion detection algorithm

disp(['Detecting motion starts and stops based on channel(s): [' num2str(idx_detect) '], with label(s): ' strjoin({EEG_motion.chanlocs(idx_detect).labels},', ')])
disp(' ')
disp('--------------------------- PLEASE CITE ----------------------------')
disp(' ')
disp('Gramann, K., Hohlefeld, F.U., Gehrke, L. et al. Human cortical dynamics during full-body heading changes.')
disp('Sci Rep 11, 18186 (2021). https://doi.org/10.1038/s41598-021-97749-8')
disp(' ')
disp('--------------------------- PLEASE CITE ----------------------------')
disp(' ')


durations = [];
latency_starts = [];
latency_stops = [];
data = EEG_motion.data(idx_detect,:);
data = sqrt(sum(data.^2,1)); % pythagoras / abs

numberstarts = 0;
numberstops = 0;

coarseThreshold = quantile(abs(data(:,search_timerange(1):search_timerange(2))),movement_threshold);

movement = false;

timePoint = search_timerange(1);
laststopTimePoint = 0;

while timePoint <= search_timerange(2)-detection_buffer_use
    step = 1;
    if ~movement
        if data(timePoint) > coarseThreshold
            
            fineThreshold = min(data(max(laststopTimePoint+1,timePoint):timePoint+detection_buffer_use))+range(data(max(laststopTimePoint+1,timePoint):timePoint+detection_buffer_use))*movement_threshold_fine;
            
            fineTimePoint = timePoint;
            
            while data(fineTimePoint) > fineThreshold && fineTimePoint > laststopTimePoint + 1
                fineTimePoint = fineTimePoint - 1;
            end
            
            numberstarts = numberstarts + 1;
            latency_starts(numberstarts) = fineTimePoint;
            movement = true;
            
            step = fineTimePoint-timePoint+1;
            
            
        end
    else
        if data(timePoint) < fineThreshold
            numberstops = numberstops + 1;
            latency_stops(numberstops) = timePoint-1;
            movement = false;
            laststopTimePoint = timePoint;
            durations(end+1) = (latency_stops(numberstops) - latency_starts(numberstarts))/EEG_motion.srate;
            
            if durations(end) <= min_duration
                durations(end) = [];
                latency_starts(numberstarts) = [];
                latency_stops(numberstops) = [];
                numberstarts = numberstarts - 1;
                numberstops = numberstops -1;
            end
        end
    end
    
    timePoint = timePoint + step;
end

latency_starts = latency_starts(1:length(durations));

disp(['Found ' num2str(length(latency_starts)) ' movements!'])

%% plot

if createplots
    disp('Plotting detection data and events')
    clear ax
    plothandles(1) = figure('color','w','position',[50 100 1500 1000]);
    ax(1) = subplot(length(idx_detect)+1,1,1);
    hold on
    plot(EEG_motion.times/1000,data)
    xlim([EEG_motion.times(1) EEG_motion.times(end)]/1000)
    plot(xlim,[coarseThreshold coarseThreshold],'k')
    title({['Detected starts (green) and stops (red) of "' eventlabel '", N = ' num2str(length(durations))]
        ['movement threshold = ' num2str(movement_threshold) ', fine movement threshold = ' num2str(movement_threshold_fine)...
        ', detection buffer = ' num2str(detection_buffer) 's, min duration = ' num2str(min_duration) 's']})
    EEG_motion.etc.motiondetect.idx_detect = idx_detect;
    EEG_motion.etc.motiondetect.eventlabel = eventlabel;
    EEG_motion.etc.motiondetect.movement_threshold = movement_threshold;
    EEG_motion.etc.motiondetect.movement_threshold_fine = movement_threshold_fine;
    EEG_motion.etc.motiondetect.detection_buffer = detection_buffer_use;
    EEG_motion.etc.motiondetect.min_duration = min_duration;
    xlabel('seconds')
    ylabel('detection data')
    
    plot(EEG_motion.times([latency_starts;latency_starts])/1000,repmat(ylim',1,length(latency_starts)),'g')
    plot(EEG_motion.times([latency_stops;latency_stops])/1000,repmat(ylim',1,length(latency_stops)),'r')
    drawnow
    
    for i_detect=1:length(idx_detect)
        disp('Plotting raw data and events')
        ax(end+1) = subplot(length(idx_detect)+1,1,1+i_detect);
        hold on
        plot(EEG_motion.times/1000,EEG_motion.data(idx_detect(i_detect),:))
        xlim([EEG_motion.times(1) EEG_motion.times(end)]/1000)
        
        title(['Channel ' num2str(idx_detect(i_detect)) ': ' EEG_motion.chanlocs(idx_detect(i_detect)).labels],'interpreter','none')
        xlabel('seconds')
        ylabel('channel data')
        
        plot(EEG_motion.times([latency_starts;latency_starts])/1000,repmat(ylim',1,length(latency_starts)),'g')
        plot(EEG_motion.times([latency_stops;latency_stops])/1000,repmat(ylim',1,length(latency_stops)),'r')
        drawnow
    end
    linkaxes(ax)
    
    plothandles(2) = figure('color','w');
    histogram(durations,0:0.05:max(durations))
    xlabel('durations (seconds)')
    ylabel('count')
    title(['Histogram of detected "' eventlabel '" durations, N = ' num2str(length(durations))])
    
end

%% add events in EEG_motion set

disp(['Adding "' eventlabel ':start" and "'  eventlabel ':stop" events to data set.'])

for i_start = 1:length(latency_starts)
    
    i_event = numel(EEG_motion.event) + 1;
    EEG_motion.event(i_event).type = [eventlabel ':start'];
    EEG_motion.event(i_event).latency = latency_starts(i_start);
    EEG_motion.event(i_event).duration = 1/EEG_motion.srate;
    
end
for i_stop = 1:length(latency_stops)
    
    i_event = numel(EEG_motion.event) + 1;
    EEG_motion.event(i_event).type = [eventlabel ':stop'];
    EEG_motion.event(i_event).latency = latency_stops(i_stop);
    EEG_motion.event(i_event).duration = 1/EEG_motion.srate;
    
end

EEG_motion = eeg_checkset(EEG_motion, 'eventconsistency');

% store info
EEG_motion.etc.motiondetect.idx_detect = idx_detect;
EEG_motion.etc.motiondetect.eventlabel = eventlabel;
EEG_motion.etc.motiondetect.latency_starts = latency_starts;
EEG_motion.etc.motiondetect.latency_stops = latency_stops;
EEG_motion.etc.motiondetect.durations = durations;
EEG_motion.etc.motiondetect.movement_threshold = movement_threshold;
EEG_motion.etc.motiondetect.movement_threshold_fine = movement_threshold_fine;
EEG_motion.etc.motiondetect.detection_buffer = detection_buffer;
EEG_motion.etc.motiondetect.min_duration = min_duration;

disp('Motion start and stop detection done!')

