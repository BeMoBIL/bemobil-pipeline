% bemobil_detect_blinks_from_ICA detects blinks and saccades based on the ICA decomposition. HEOG and VEOG ICs are automatically
% detected based on their topographies and their spectral power <5Hz. Blinks and saccades are detected using findpeaks
% with distance, prominence, and width of the peaks. The parameters can be set and an informed decision about the
% detector efficacy can be made by the plots of the detection which includes a histogram of the prominences and widths
% (including those exceeding the thresholds). 
%
% All parameters, the found blink and saccade latencies, as well as the prominences and widths of the detector are
% stored in the EEG.etc.ICs4events.eye field.
%
% The algorithm was first used and published in 
%
% Input:
%       EEG                             - EEGLAB dataset containing an ICA decomposition
%       idx_VEOG                        - OPTIONAL index of the VEOG IC (default = will be identified automatically)
%       idx_HEOG                        - OPTIONAL index of the HEOG IC (default = will be identified automatically)
%       medianFilterLength              - OPTIONAL filter duration of the median filter for the eye data in seconds (default = 0.08)
%       blink_minPeakDistance           - OPTIONAL findpeaks min peak distance in seconds (default = 0.15)
%       blink_minPeakWidth              - OPTIONAL findpeaks min peak width in seconds (default = 0.15)
%       blink_maxPeakWidth              - OPTIONAL findpeaks min peak distance in seconds (default = 0.3)
%       blink_minProminence             - OPTIONAL findpeaks min peak promincences in uV (default = 20)
%       saccade_minPeakDistance         - OPTIONAL findpeaks min peak distance in seconds (default = 0.15)
%       saccade_minPeakWidth            - OPTIONAL findpeaks min peak width in seconds (default = 0)
%       saccade_maxPeakWidth            - OPTIONAL findpeaks min peak distance in seconds (default = 0.3)
%       saccade_minProminence           - OPTIONAL findpeaks min peak promincences in uV^2 (default = 4)
%       saccade_blink_distance          - OPTIONAL distance a saccade needs to have to a blink, all closer will be
%                                           deleted, in seconds (default = 0.1)
%       store_saccades                  - OPTIONAL whether or not saccades should be stored as events in the end (default = 1)
%
% Output:
%       EEG                 - EEGLAB dataset containing detected events
%       plothandles         - handle to the created plots to allow script-based saving and closing
%
% Example:
%        [EEG, plothandles] = bemobil_detect_blinks_from_ICA(EEG)
%        [EEG, plothandles] = bemobil_detect_blinks_from_ICA(EEG,'blink_minProminence',30)
%
% See also: findpeaks
%
% Authors: Marius Klug, Anna Wunderlich, 2022

function [EEG, plothandles] = bemobil_detect_blinks_from_ICA(EEG, varargin)

if nargin == 0
    help bemobil_detect_blinks_from_ICA
    return
end

p = inputParser;

% set the desired and optional input arguments
addRequired(p, 'EEG', @(x) validateattributes(x,{'struct'},{},'bemobil_detect_blinks_from_ICA','EEG'));
addOptional(p, 'idx_VEOG', [],  @(x) validateattributes(x,{'numeric'},{'positive','integer','scalar'},'bemobil_detect_blinks_from_ICA','idx_VEOG')) 
addOptional(p, 'idx_HEOG', [],  @(x) validateattributes(x,{'numeric'},{'positive','integer','scalar'},'bemobil_detect_blinks_from_ICA','idx_HEOG')) 
addOptional(p, 'medianFilterLength', 0.08,  @(x) validateattributes(x,{'numeric'},{'positive','scalar'},'bemobil_detect_blinks_from_ICA','medianFilterLength')) 
addOptional(p, 'blink_minPeakDistance', 0.15,  @(x) validateattributes(x,{'numeric'},{'positive','scalar'},'bemobil_detect_blinks_from_ICA','blink_minPeakDistance')) 
addOptional(p, 'blink_minPeakWidth', 0.04,  @(x) validateattributes(x,{'numeric'},{'positive','scalar'},'bemobil_detect_blinks_from_ICA','blink_minPeakWidth')) 
addOptional(p, 'blink_maxPeakWidth', 0.3,  @(x) validateattributes(x,{'numeric'},{'positive','scalar'},'bemobil_detect_blinks_from_ICA','blink_maxPeakWidth')) 
addOptional(p, 'blink_minProminence', 20,  @(x) validateattributes(x,{'numeric'},{'positive','scalar'},'bemobil_detect_blinks_from_ICA','blink_minProminence')) 
addOptional(p, 'saccade_minPeakDistance', 0.2,  @(x) validateattributes(x,{'numeric'},{'positive','scalar'},'bemobil_detect_blinks_from_ICA','saccade_minPeakDistance')) 
addOptional(p, 'saccade_minPeakWidth', 0,  @(x) validateattributes(x,{'numeric'},{'positive','scalar'},'bemobil_detect_blinks_from_ICA','saccade_minPeakWidth')) 
addOptional(p, 'saccade_minProminence', 4,  @(x) validateattributes(x,{'numeric'},{'positive','scalar'},'bemobil_detect_blinks_from_ICA','saccade_minProminence')) 
addOptional(p, 'saccade_blink_distance', 0.1,  @(x) validateattributes(x,{'numeric'},{'positive','scalar'},'bemobil_detect_blinks_from_ICA','saccade_blink_distance')) 
addOptional(p, 'store_saccades', 1,  @(x) validateattributes(x,{'numeric','logical'},{'scalar','binary'},'bemobil_detect_blinks_from_ICA','store_saccades')) 


% parse the input
parse(p,EEG,varargin{:});

% then set/get all the inputs out of this structure
EEG = p.Results.EEG;
idx_VEOG = p.Results.idx_VEOG;
idx_HEOG = p.Results.idx_HEOG;
medianFilterLength = p.Results.medianFilterLength;
blink_minPeakDistance = p.Results.blink_minPeakDistance;
blink_minPeakWidth = p.Results.blink_minPeakWidth;
blink_maxPeakWidth = p.Results.blink_maxPeakWidth;
blink_minProminence = p.Results.blink_minProminence;
saccade_minPeakDistance = p.Results.saccade_minPeakDistance;
saccade_minPeakWidth = p.Results.saccade_minPeakWidth;
saccade_minProminence = p.Results.saccade_minProminence;
saccade_blink_distance = p.Results.saccade_blink_distance;
store_saccades = p.Results.store_saccades;

if ~isfield(EEG,'data')
    error('Wrong data type for EEG input! Needs to be an EEGLAB struct.')
end

EEG.etc.ICs4events.eye = [];
EEG.etc.ICs4events.eye.idx_VEOG = idx_VEOG;
EEG.etc.ICs4events.eye.idx_HEOG = idx_HEOG;
EEG.etc.ICs4events.eye.medianFilterLength = medianFilterLength;
EEG.etc.ICs4events.eye.blink_minPeakDistance = blink_minPeakDistance;
EEG.etc.ICs4events.eye.blink_minPeakWidth = blink_minPeakWidth;
EEG.etc.ICs4events.eye.blink_maxPeakWidth = blink_maxPeakWidth;
EEG.etc.ICs4events.eye.blink_minProminence = blink_minProminence;
EEG.etc.ICs4events.eye.saccade_minPeakDistance = saccade_minPeakDistance;
EEG.etc.ICs4events.eye.saccade_minPeakWidth = saccade_minPeakWidth;
EEG.etc.ICs4events.eye.saccade_minProminence = saccade_minProminence;
EEG.etc.ICs4events.eye.saccade_blink_distance = saccade_blink_distance;

%% search for EOG indices

if ~isfield(EEG.etc,'ic_classification')
    disp('No ICLabel info found. Classifying ICs using ICLabel defaults.')
    EEG = iclabel(EEG);
end

% search HEOG and VEOG ICs using ICLabel
eye_ICs = find(max(EEG.etc.ic_classification.ICLabel.classifications,[],2)==...
    EEG.etc.ic_classification.ICLabel.classifications(:,3))';

idx_VEOG = 0;
invert_veog = 0;
power_veog = 0;

idx_HEOG = 0;
sum_heog = 0;

plothandles(1) = figure('position',[285 202 790 420]);

plot_idx = 0;
for idx = eye_ICs
    
    % plot
    plot_idx = plot_idx+1;
    subplot(1,length(eye_ICs),plot_idx)
    [~, topoimage]= topoplot(EEG.icawinv(:,idx), EEG.chanlocs,'noplot','off');
    
    % simple classifier for heog and veog because they are really easy to distinguish
    this_sum_veog = sum(sum(topoimage(50:60,30:40)));
    this_sum_heog = sum(sum(topoimage(48:52,15:20))) - sum(sum(topoimage(48:52,50:55)));
    
    % check low spectral power
    [pxx, f] = pwelch(EEG.icaact(idx,:),EEG.srate*10,[],[],EEG.srate); 
    low_power = mean(pxx(f<5));
        
    title({['IC #' num2str(num2str(idx))]
        ['sum HEOG: ' num2str(abs(this_sum_heog))]
        ['sum VEOG: ' num2str(abs(this_sum_veog))]
        ['power <5Hz: ' num2str(low_power)]})
    
    % check corr with HEOG prototype
    if abs(this_sum_heog)>abs(this_sum_veog) && abs(this_sum_heog) > sum_heog
        idx_HEOG = idx;
        sum_heog = abs(this_sum_heog);
    elseif low_power > power_veog
        idx_VEOG = idx;
        power_veog = low_power;
        if this_sum_veog < 0
            invert_veog = 1;
        end
    end
end

% save selected ICs in EEG.etc
if isempty(EEG.etc.ICs4events.eye.idx_HEOG)
    EEG.etc.ICs4events.eye.idx_HEOG = idx_HEOG;
else
    disp('HEOG was provided, found HEOG will be discarded.');
end
if isempty(EEG.etc.ICs4events.eye.idx_VEOG)
    EEG.etc.ICs4events.eye.idx_VEOG = idx_VEOG;
else
    disp('HEOG was provided, found HEOG will be discarded.');
end

%% detect blinks

temp_VEOG = EEG.icaact(EEG.etc.ICs4events.eye.idx_VEOG,:);
temp_HEOG = EEG.icaact(EEG.etc.ICs4events.eye.idx_HEOG,:);

% median filter
smoothVEOG = smoothdata(temp_VEOG,'movmedian',round(medianFilterLength*EEG.srate));
smoothHEOG = smoothdata(temp_HEOG,'movmedian',round(medianFilterLength*EEG.srate));

if invert_veog
    smoothVEOG = -smoothVEOG;
end

% initialize and define all parameters for blink detection
minPeakDistance = blink_minPeakDistance;
minPeakWidth = blink_minPeakWidth;
maxPeakWidth = blink_maxPeakWidth;
minPeakProminence = blink_minProminence;

% now find the actual blinks without adding the time so the locs are correct
[~,blinkLatencies] = findpeaks(smoothVEOG,'MinPeakProminence',minPeakProminence,...
    'MinPeakDistance',minPeakDistance*EEG.srate,...
    'MinPeakWidth',minPeakWidth*EEG.srate,'MaxPeakWidth',maxPeakWidth*EEG.srate);

% find again without limits for width and prominence for plot, using EEG.times to be easier interpretable
[~,~,w,p] = findpeaks(smoothVEOG,EEG.times/1000,'MinPeakDistance',minPeakDistance);

% find again for plot only, using EEG.times to be easier interpretable
plothandles(2) = figure('position',[111 365 1083 419],'color','w');
subplot(2,2,1:2);
findpeaks(smoothVEOG,EEG.times/1000,'MinPeakProminence',minPeakProminence,...
    'MinPeakDistance',minPeakDistance,'annotate','extents',...
    'MinPeakWidth',minPeakWidth,'MaxPeakWidth',maxPeakWidth)
title(['blinks, IC# ' num2str(EEG.etc.ICs4events.eye.idx_VEOG)])
xlabel('time [s]')
ylabel('smoothed VEOG [\muV]')

% plot prominences and widths
subplot(223);
histogram(p,'binwidth',0.5)
title('all prominences')
hold on
l=plot([minPeakProminence minPeakProminence],ylim,'r');
legend(l,'min prominence')
xlabel('\muV')

subplot(224);
histogram(w,'binwidth',0.005)
title('all widths')
hold on
l(1)=plot([minPeakWidth minPeakWidth],ylim,'r');
l(2)=plot([maxPeakWidth maxPeakWidth],ylim,'k');
legend(l,{'min width','max width'})
xlabel('seconds')

% store
EEG.etc.ICs4events.eye.blink_prominences = p;
EEG.etc.ICs4events.eye.blink_widths = w;

%% detect saccades

% compute REOG
REOG = sqrt(smoothVEOG.^2 + smoothHEOG.^2);

% compute the squared derivatie
sqdiff_REOG = diff(REOG).^2;

% initialize all variables and parameters necessary for saccade events
minPeakDistance = saccade_minPeakDistance;
minPeakWidth = saccade_minPeakWidth;
maxPeakWidth = 9999999;
minPeakProminence = saccade_minProminence;

% find peaks in the squared derivative of REOG
[~,sacLatencies] = findpeaks(sqdiff_REOG,'MinPeakProminence',minPeakProminence,...
    'MinPeakDistance',minPeakDistance*EEG.srate,...
    'MinPeakWidth',minPeakWidth*EEG.srate,'MaxPeakWidth',maxPeakWidth*EEG.srate);

% find again without limits for width and prominence for plot
[~,~,w,p] = findpeaks(sqdiff_REOG,EEG.times(2:end)/1000,'MinPeakDistance',minPeakDistance);

% find again for plot only
plothandles(3) = figure('position',[111 365 1083 419],'color','w');
subplot(2,2,1:2);
findpeaks(sqdiff_REOG,EEG.times(2:end)/1000,'MinPeakProminence',minPeakProminence,...
    'MinPeakDistance',minPeakDistance,'annotate','extents',...
    'MinPeakWidth',minPeakWidth,'MaxPeakWidth',maxPeakWidth);
title(['saccades, IC #' num2str(EEG.etc.ICs4events.eye.idx_HEOG)])
xlabel('time [s]')
ylabel({'squared derivative of','combined EOGs [\muV^2]'})

% plot prominences and widths
subplot(223);
histogram(p,'binwidth',0.5)
title('all prominences')
hold on
l=plot([minPeakProminence minPeakProminence],ylim,'r');
legend(l,'min prominence')
xlabel('\muV^2')

subplot(224);
histogram(w,'binwidth',0.0005)
title('all widths')
xlabel('seconds')
hold on
l(1)=plot([minPeakWidth minPeakWidth],ylim,'r');
legend(l,{'min width'})

% store
EEG.etc.ICs4events.eye.sacccade_prominences = p;
EEG.etc.ICs4events.eye.sacccade_widths = w;

% check whether a detected saccade is too close to a blink and likely a blink
deleteSac = [];
for sacIndex = 1:length(sacLatencies)
    if min(abs(blinkLatencies - sacLatencies(sacIndex))) < saccade_blink_distance*EEG.srate
        deleteSac(end+1) = sacIndex;
    end
end

% delete those saccades
sacLatencies(deleteSac) = [];

%% add events
EEG.etc.ICs4events.eye.blinkLatencies = blinkLatencies;
EEG.etc.ICs4events.eye.sacLatencies = sacLatencies;

% create blink events
event_latencies = blinkLatencies;
for latency = event_latencies
    i = numel(EEG.event) + 1;
    EEG.event(i).type = 'blink';
    EEG.event(i).latency = latency;
    EEG.event(i).duration = 1/EEG.srate;
end

if store_saccades
    % create saccade events
    event_latencies = sacLatencies;
    for latency = event_latencies
        i = numel(EEG.event) + 1;
        EEG.event(i).type = 'saccade';
        EEG.event(i).latency = latency;
        EEG.event(i).duration = 1/EEG.srate;
    end
end

EEG = eeg_checkset(EEG, 'eventconsistency');

%% plot final eye data only

rejcomps = 1:size(EEG.icaact,1);
rejcomps([EEG.etc.ICs4events.eye.idx_VEOG EEG.etc.ICs4events.eye.idx_HEOG]) = [];

EEG_plot = pop_subcomp( EEG, rejcomps, 0);
pop_eegplot( EEG_plot, 0, 1, 1); 
plothandles(end+1) = gcf;

disp('Blink and saccade detection done!')

