% bemobil_gait_analysis detects step cycle events from motion capture of the feet. Detected events for both feet are:
% 
% foot movement start (= flatfoot phase end) -> toe-off -> heelstrike -> foot movement stop (= flatfoot phase start)
% 
% Foot movements are detected using the function bemobil_detect_motion_startstops on the z-Axis of the tracking (up-down
% movement). Toe-off and heelstrike events are then defined by the zero-crossing of the velocity in the x-Axis
% (forward-backward), after, or before the foot movement events, respectively. For this both x and y axes need to be
% provided and an artificial coordinate system is constructed using a PCA, since the coordinate system is not always
% directed at the front.
%
% The tracking needs to be either on a treadmill (stationary movement), or relative to an upper body part like the torso
% or head, in order to provide foot movement that goes forward-backward in each step. Hence, if the participant did walk
% on a treadmill, only feet tracking channel indices are necessary, if the experiment contained natural overground
% walking, the x and y indices of the head (or another upper torso tracker) are necessary to obtain correct events.
%
% All parameters and the found step event latencies are stored in the EEG_motion.etc.gait_analysis field.
%
% Input:
%       EEG_motion                          - EEGLAB dataset containing motion capture
%       idx_x_1                             - index of the first foot x-Axis (forward-backward) (REQUIRED)
%       idx_y_1                             - index of the first foot y-Axis (left-right) (REQUIRED)
%       idx_z_1                             - index of the first foot z-Axis (up-down) (REQUIRED)
%       idx_x_2                             - index of the second foot x-Axis (forward-backward) (REQUIRED)
%       idx_y_2                             - index of the second foot y-Axis (left-right) (REQUIRED)
%       idx_z_2                             - index of the second foot z-Axis (up-down) (REQUIRED)
%       search_timerange                    - OPTIONAL latency timerange where the gait should be analyzed (in samples, default = all)
%       idx_head_x                          - OPTIONAL index of the head x-Axis (forward-backward)
%       idx_head_y                          - OPTIONAL index of the head y-Axis (left-right)
%       movement_threshold_quantile         - OPTIONAL movement detection threshold (default = 0.5)
%       movement_threshold_fine_quantile    - OPTIONAL fine movement detection threshold (default = 0.05)
%       min_duration                        - OPTIONAL minimum foot movement duration (in s, default = 0.5)
%       detection_buffer                    - OPTIONAL fine movement detection buffer (in s, default = 1)
%
% Output:
%       EEG_motion                          - EEGLAB dataset containing detected events
%       plothandles                         - handle to the created plots to allow script-based saving and closing
%
% Example:
%        [EEG_motion, plothandles] = bemobil_detect_blinks(EEG_motion, idx_x_1, idx_y_1, idx_z_1, idx_x_2, idx_y_2, idx_z_2)
%        [EEG_motion, plothandles] = bemobil_detect_blinks(EEG_motion, idx_x_1, idx_y_1, idx_z_1, idx_x_2, idx_y_2, idx_z_2,...
%                                       'search_timerange',search_timerange,'idx_head_x',idx_head_x,'idx_head_y',idx_head_y)
%
% See also: bemobil_detect_motion_startstops
%
% Authors: Marius Klug, 2022

function [EEG_motion, plothandles] = bemobil_gait_analysis(EEG_motion, idx_x_1, idx_y_1, idx_z_1, idx_x_2, idx_y_2, idx_z_2, varargin)

if nargin == 0
    help bemobil_gait_analysis
    return
end

p = inputParser;

% set the desired and optional input arguments
addRequired(p, 'EEG_motion', @(x) validateattributes(x,{'struct'},{},'bemobil_gait_analysis','EEG_motion'));
addRequired(p, 'idx_x_1',  @(x) validateattributes(x,{'numeric'},{'positive','integer','scalar'},'bemobil_gait_analysis','idx_x_1'))
addRequired(p, 'idx_y_1',  @(x) validateattributes(x,{'numeric'},{'positive','integer','scalar'},'bemobil_gait_analysis','idx_y_1'))
addRequired(p, 'idx_z_1',  @(x) validateattributes(x,{'numeric'},{'positive','integer','scalar'},'bemobil_gait_analysis','idx_z_1'))
addRequired(p, 'idx_x_2',  @(x) validateattributes(x,{'numeric'},{'positive','integer','scalar'},'bemobil_gait_analysis','idx_x_2'))
addRequired(p, 'idx_y_2',  @(x) validateattributes(x,{'numeric'},{'positive','integer','scalar'},'bemobil_gait_analysis','idx_y_2'))
addRequired(p, 'idx_z_2',  @(x) validateattributes(x,{'numeric'},{'positive','integer','scalar'},'bemobil_gait_analysis','idx_z_2'))
addOptional(p, 'search_timerange', [],  @(x) validateattributes(x,{'numeric'},{'positive','vector'},'bemobil_gait_analysis','search_timerange'))
addOptional(p, 'idx_head_x', [],  @(x) validateattributes(x,{'numeric'},{'positive','integer','scalar'},'bemobil_gait_analysis','idx_head_x'))
addOptional(p, 'idx_head_y', [],  @(x) validateattributes(x,{'numeric'},{'positive','integer','scalar'},'bemobil_gait_analysis','idx_head_x'))
addOptional(p, 'movement_threshold_quantile', 0.5,  @(x) validateattributes(x,{'numeric'},{'positive','scalar'},'bemobil_gait_analysis','movement_threshold_quantile'))
addOptional(p, 'movement_threshold_fine_quantile', 0.05,  @(x) validateattributes(x,{'numeric'},{'positive','scalar'},'bemobil_gait_analysis','movement_threshold_fine_quantile'))
addOptional(p, 'min_duration', 0.5,  @(x) validateattributes(x,{'numeric'},{'positive','scalar'},'bemobil_gait_analysis','min_duration'))
addOptional(p, 'detection_buffer', 1,  @(x) validateattributes(x,{'numeric'},{'positive','scalar'},'bemobil_gait_analysis','detection_buffer'))

% parse the input
parse(p,EEG_motion,idx_x_1,idx_y_1,idx_z_1,idx_x_2,idx_y_2,idx_z_2,varargin{:});

% then set/get all the inputs out of this structure
EEG_motion = p.Results.EEG_motion;
idx_x_1 = p.Results.idx_x_1;
idx_y_1 = p.Results.idx_y_1;
idx_z_1 = p.Results.idx_z_1;
idx_x_2 = p.Results.idx_x_2;
idx_y_2 = p.Results.idx_y_2;
idx_z_2 = p.Results.idx_z_2;
search_timerange = p.Results.search_timerange;
idx_head_x = p.Results.idx_head_x;
idx_head_y = p.Results.idx_head_y;
movement_threshold_quantile = p.Results.movement_threshold_quantile;
movement_threshold_fine_quantile = p.Results.movement_threshold_fine_quantile;
min_duration = p.Results.min_duration;
detection_buffer = p.Results.detection_buffer;

if ~isfield(EEG_motion,'data')
    error('Wrong data type for EEG_motion input! Needs to be an EEGLAB struct.')
end

if isempty(search_timerange)
    search_timerange = [1 EEG_motion.pnts];
end

if xor(isempty(idx_head_x),isempty(idx_head_y))
    error('If head indices are used, both x and y must be set.')
end

% allow several gait detectors in succession
if isfield(EEG_motion.etc,'gait_analysis')
    idx_detection = length(EEG_motion.etc.gait_analysis) + 1;
else
    idx_detection = 1;
end

EEG_motion.etc.gait_analysis(idx_detection).idx_x_1 = idx_x_1;
EEG_motion.etc.gait_analysis(idx_detection).idx_y_1 = idx_y_1;
EEG_motion.etc.gait_analysis(idx_detection).idx_z_1 = idx_z_1;
EEG_motion.etc.gait_analysis(idx_detection).idx_x_2 = idx_x_2;
EEG_motion.etc.gait_analysis(idx_detection).idx_y_2 = idx_y_2;
EEG_motion.etc.gait_analysis(idx_detection).idx_z_2 = idx_z_2;
EEG_motion.etc.gait_analysis(idx_detection).search_timerange = search_timerange;
EEG_motion.etc.gait_analysis(idx_detection).idx_head_x = idx_head_x;
EEG_motion.etc.gait_analysis(idx_detection).idx_head_y = idx_head_y;
EEG_motion.etc.gait_analysis(idx_detection).movement_threshold_quantile = movement_threshold_quantile;
EEG_motion.etc.gait_analysis(idx_detection).movement_threshold_fine_quantile = movement_threshold_fine_quantile;
EEG_motion.etc.gait_analysis(idx_detection).min_duration = min_duration;
EEG_motion.etc.gait_analysis(idx_detection).detection_buffer = detection_buffer;


%% detect general foot movements

% lowpass for z axis movement detection only
EEG_motion_filt = pop_eegfiltnew(EEG_motion, 0.5, [], [], 0, [], 0);
% make all data positive for the detector to work
EEG_motion_filt.data(idx_z_1,:) = EEG_motion_filt.data(idx_z_1,:)-min(EEG_motion_filt.data(idx_z_1,:));
EEG_motion_filt.data(idx_z_2,:) = EEG_motion_filt.data(idx_z_2,:)-min(EEG_motion_filt.data(idx_z_2,:));

[EEG_motion_filt, plots1] = bemobil_detect_motion_startstops(EEG_motion_filt,idx_z_1,'foot1',movement_threshold_quantile,movement_threshold_fine_quantile,min_duration,detection_buffer,search_timerange,1);
[EEG_motion_filt, plots2] = bemobil_detect_motion_startstops(EEG_motion_filt,idx_z_2,'foot2',movement_threshold_quantile,movement_threshold_fine_quantile,min_duration,detection_buffer,search_timerange,1);

EEG_motion.event = EEG_motion_filt.event;
plothandles = [plots1 plots2];

%% create data x y data using PCA for both feet 
clear scores
for i_foot = 1:2
    
    if i_foot==1
        idx_x = idx_x_1;
        idx_y = idx_y_1;
    else
        idx_x = idx_x_2;
        idx_y = idx_y_2;
    end
    
    % make data head-relative if selected
    if isempty(idx_head_x)
        xy_data = EEG_motion.data([idx_x idx_y],:);
    else
        xy_data = EEG_motion.data([idx_x idx_y],:)- EEG_motion.data([idx_head_x idx_head_y],:);
    end
    
    
    % PCA of the x and y axis to find the real gait parameters
    [coeff,~,~,~,~,mu] = pca(xy_data(:,search_timerange(1):search_timerange(2))');
    
    scores(i_foot,1:2,1:length(xy_data)) = transpose((xy_data'-mu)*coeff);
end

%% find toeoff and heelstrike events 

gaitspeedmediansmoothing = 25;

% get forward-backward positions and velocity
pca_x = squeeze(scores(:,1,:));
pca_x_vel = diff(pca_x,1,2);

latencies_toeoff = NaN(2,99999);
latencies_heelstrike = NaN(2,99999);
latencies_footstart = NaN(2,99999);
latencies_footend = NaN(2,99999);

% !!!!!!!!!!!!!!!!!!!!!
% these gait parameters are experiemntal and currently not in use
% !!!!!!!!!!!!!!!!!!!!!
gaitspeeds = NaN(2,99999);
smoothed_gaitspeeds = NaN(2,99999);
stride_lengths = NaN(2,99999);

i_latencies_toeoff = [0 0];
i_latencies_heelstrike = [0 0];
i_latencies_footstart = [0 0];
i_latencies_footend = [0 0];
i_gaitspeeds = [0 0];
i_smoothed_gaitspeeds = [0 0];
i_stride_lengths = [0 0];

for i=1:length(EEG_motion.event) % go through all events
    
    for i_foot = 1:2
        
        eventstring = ['foot' num2str(i_foot)];
        if strcmp(EEG_motion.event(i).type,[eventstring ':start']) % foot movement start (= flat foot end)
            i_latencies_footstart(i_foot) = i_latencies_footstart(i_foot)+1;
            latencies_footstart(i_foot,i_latencies_footstart(i_foot)) = EEG_motion.event(i).latency;
            
            % from foot movement start go forward to find the zero crossing, which is taken as toe-off
            % IMPORTANT: for this to work the foot x/y must be stationary (treadmill or relative to head)
            i_toeoff = EEG_motion.event(i).latency;
            while sign(pca_x_vel(i_foot,i_toeoff)) ==...
                    sign(pca_x_vel(i_foot,latencies_footstart(i_foot,i_latencies_footstart(i_foot))))
                i_toeoff = i_toeoff + 1;
            end
            i_latencies_toeoff(i_foot) = i_latencies_toeoff(i_foot)+1;
            latencies_toeoff(i_foot,i_latencies_toeoff(i_foot)) = i_toeoff;
            
            % store gait speed during previous flat foot phase (sum of velocity divided by duration) in m/s
            if i_latencies_footend(i_foot)>0
                i_gaitspeeds(i_foot) = i_gaitspeeds(i_foot)+1;
                gaitspeeds(i_foot,i_gaitspeeds(i_foot)) =...
                    abs(sum(pca_x_vel(i_foot,latencies_footend(i_foot,i_latencies_footend(i_foot)):latencies_footstart(i_foot,i_latencies_footstart(i_foot)))) /...
                    ((latencies_footstart(i_foot,i_latencies_footstart(i_foot))-latencies_footend(i_foot,i_latencies_footend(i_foot)))/EEG_motion.srate));
            end
            
        elseif strcmp(EEG_motion.event(i).type,[eventstring ':stop']) % foot movement stop (= flat foot start)
            i_latencies_footend(i_foot) = i_latencies_footend(i_foot)+1;
            latencies_footend(i_foot,i_latencies_footend(i_foot)) = EEG_motion.event(i).latency;
            
            % from foot movement stop go backward to find the zero crossing, which is taken as heelstrike
            i_heelstrike = EEG_motion.event(i).latency;
            while sign(pca_x_vel(i_foot,i_heelstrike)) ==...
                    sign(pca_x_vel(i_foot,latencies_footend(i_foot,i_latencies_footend(i_foot))))
                i_heelstrike = i_heelstrike - 1;
            end
            i_latencies_heelstrike(i_foot) = i_latencies_heelstrike(i_foot)+1;
            latencies_heelstrike(i_foot,i_latencies_heelstrike(i_foot)) = i_heelstrike;
            
            % store gait parameters
            if i_gaitspeeds(i_foot)>0 && i_latencies_heelstrike(i_foot)>1
                
                % reduce variability of gait speed a bit to be more consistent, especially get rid of outliers due to
                % measurement issues
                i_smoothed_gaitspeeds(i_foot) = i_smoothed_gaitspeeds(i_foot)+1;
                smoothed_gaitspeeds(i_foot,i_smoothed_gaitspeeds(i_foot)) =...
                    median(gaitspeeds(i_foot,max(i_gaitspeeds(i_foot)-gaitspeedmediansmoothing+1,1):i_gaitspeeds(i_foot)));
                
                % stride length is the heelstrike to heelstrike distance, which is measured in the x difference plus the
                % traveled distance according to the velocity in the meantime (in m)
                i_stride_lengths(i_foot) = i_stride_lengths(i_foot)+1;
                stride_lengths(i_foot,i_stride_lengths(i_foot)) =...
                    (pca_x(i_foot,latencies_heelstrike(i_foot,i_latencies_heelstrike(i_foot))) -...
                    pca_x(i_foot,latencies_heelstrike(i_foot,i_latencies_heelstrike(i_foot)-1)))*...
                    sign(pca_x(i_foot,latencies_heelstrike(i_foot,i_latencies_heelstrike(i_foot)-1))) + ... % pca is agnostic to forward/backward, so we define this to be invariant
                    ((latencies_heelstrike(i_foot,i_latencies_heelstrike(i_foot)) -...
                    latencies_heelstrike(i_foot,i_latencies_heelstrike(i_foot)-1))/EEG_motion.srate * smoothed_gaitspeeds(end));
                
                % step length is the distance one foot travels before the other (heelstrike to heelstrike)
%                 step_lengths(end+1)
            end
        end
    end
end

latencies_toeoff(:,all(isnan(latencies_toeoff))) = [];
latencies_heelstrike(:,all(isnan(latencies_heelstrike))) = [];
latencies_footstart(:,all(isnan(latencies_footstart))) = [];
latencies_footend(:,all(isnan(latencies_footend))) = [];
gaitspeeds(:,all(isnan(gaitspeeds))) = [];
smoothed_gaitspeeds(:,all(isnan(smoothed_gaitspeeds))) = [];
stride_lengths(:,all(isnan(stride_lengths))) = [];

%% store events 

EEG_motion.etc.gait_analysis.latencies_footstart = latencies_footstart;
EEG_motion.etc.gait_analysis.latencies_toeoff = latencies_toeoff;
EEG_motion.etc.gait_analysis.latencies_heelstrike = latencies_heelstrike;
EEG_motion.etc.gait_analysis.latencies_footend = latencies_footend;

for i_foot = 1:2
    eventstring = ['foot' num2str(i_foot)];
    for latency = latencies_toeoff(i_foot,~isnan(latencies_toeoff(i_foot,:)))
        i = numel(EEG_motion.event) + 1;
        EEG_motion.event(i).type = [eventstring ':toeoff'];
        EEG_motion.event(i).latency = latency;
        EEG_motion.event(i).duration = 1/EEG_motion.srate;
    end
    for latency = latencies_heelstrike(i_foot,~isnan(latencies_heelstrike(i_foot,:)))
        i = numel(EEG_motion.event) + 1;
        EEG_motion.event(i).type = [eventstring ':heelstrike'];
        EEG_motion.event(i).latency = latency;
        EEG_motion.event(i).duration = 1/EEG_motion.srate;
    end
end
    
EEG_motion = eeg_checkset(EEG_motion, 'eventconsistency');

%% plot final steps only
alllabels = {EEG_motion.chanlocs.labels};

EEG_motion_plot = pop_select( EEG_motion, 'channel',alllabels([idx_x_1, idx_y_1, idx_z_1, idx_x_2, idx_y_2, idx_z_2]));
pop_eegplot( EEG_motion_plot, 1, 1, 1); 
plothandles(end+1) = gcf;

disp('Gait analysis done!')
