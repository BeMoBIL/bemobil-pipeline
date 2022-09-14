% bemobil_detect_steps_from_ICA detects walking steps based on the ICA decomposition and a provided gait IC index. Steps are
% detected using findpeaks with distance, prominence, and width of the peaks. The parameters can be set and
% an informed decision about the detector efficacy can be made by the plots of the detection which includes a histogram
% of the prominences and widths (including those exceeding the thresholds). The detection can be repeated with varying
% search boundaries or IC indices.
%
% All parameters, the found step latencies, as well as the prominences and widths of the detector are
% stored in the EEG.etc.ICs4events.gait.
%
% The algorithm was first used and published in 
%
% Input:
%       EEG                             - EEGLAB dataset containing an ICA decomposition
%       idx_gait                        - index of the step IC (REQUIRED)
%       startSearch                     - OPTIONAL latency when the search should start (in samples, default = 1)
%       endSearch                       - OPTIONAL latency when the search should end (in samples, default = all)
%       medianFilterLength              - OPTIONAL filter duration of the median filter for the gait data in seconds (default = 0.1)
%       minStepDistance                 - OPTIONAL findpeaks min peak distance in seconds (default = 0.4)
%       minPeakDuration                 - OPTIONAL findpeaks min peak width in seconds (default = 0.05)
%       maxPeakDuration                 - OPTIONAL findpeaks max peak width in seconds (default = 999)
%       minPeakProminence               - OPTIONAL findpeaks min prominence in uV (default = 0)
%
% Output:
%       EEG                 - EEGLAB dataset containing detected events
%       plothandles         - handle to the created plots to allow script-based saving and closing
%
% Example:
%        [EEG, plothandles] = bemobil_detect_steps_from_ICA(EEG,idx_step_IC)
%        [EEG, plothandles] = bemobil_detect_steps_from_ICA(EEG,idx_step_IC,'endSearch',round(3373*EEG.srate))
%
% See also: findpeaks
%
% Authors: Marius Klug, Anna Wunderlich, 2022

function [EEG, plothandles] = bemobil_detect_steps_from_ICA(EEG, idx_gait, varargin)

if nargin == 0
    help bemobil_detect_steps_from_ICA
    return
end

p = inputParser;

% set the desired and optional input arguments
addRequired(p, 'EEG', @(x) validateattributes(x,{'struct'},{},'bemobil_detect_blinks','EEG'));
addRequired(p, 'idx_gait',  @(x) validateattributes(x,{'numeric'},{'positive','integer','scalar'},'bemobil_detect_blinks','idx_gait')) 
addOptional(p, 'startSearch', 1,  @(x) validateattributes(x,{'numeric'},{'positive','scalar'},'bemobil_detect_blinks','startSearch')) 
addOptional(p, 'endSearch', 0,  @(x) validateattributes(x,{'numeric'},{'positive','scalar'},'bemobil_detect_blinks','endSearch')) 
addOptional(p, 'medianFilterLength', 0.1,  @(x) validateattributes(x,{'numeric'},{'positive','scalar'},'bemobil_detect_blinks','medianFilterLength')) 
addOptional(p, 'minStepDistance', 0.4,  @(x) validateattributes(x,{'numeric'},{'positive','scalar'},'bemobil_detect_blinks','minStepDistance')) 
addOptional(p, 'minPeakDuration', 0.05,  @(x) validateattributes(x,{'numeric'},{'positive','scalar'},'bemobil_detect_blinks','minPeakDuration')) 
addOptional(p, 'maxPeakDuration', 999,  @(x) validateattributes(x,{'numeric'},{'positive','scalar'},'bemobil_detect_blinks','maxPeakDuration')) 
addOptional(p, 'minPeakProminence', 0,  @(x) validateattributes(x,{'numeric'},{'positive','scalar'},'bemobil_detect_blinks','minPeakProminence')) 

% parse the input
parse(p,EEG,idx_gait,varargin{:});

% then set/get all the inputs out of this structure
EEG = p.Results.EEG;
idx_gait = p.Results.idx_gait;
startSearch = p.Results.startSearch;
endSearch = p.Results.endSearch;
medianFilterLength = p.Results.medianFilterLength;
minStepDistance = p.Results.minStepDistance;
minPeakDuration = p.Results.minPeakDuration;
maxPeakDuration = p.Results.maxPeakDuration;
minPeakProminence = p.Results.minPeakProminence;

if ~isfield(EEG,'data')
    error('Wrong data type for EEG input! Needs to be an EEGLAB struct.')
end

if endSearch == 0
    endSearch = EEG.pnts;
end

% allow several gait detectors in succession
if isfield(EEG.etc,'ICs4events') && isfield(EEG.etc.ICs4events,'gait')
    idx_detection = length(EEG.etc.ICs4events.gait) + 1;
else
    idx_detection = 1;
end

EEG.etc.ICs4events.gait(idx_detection).idx_gait = idx_gait;
EEG.etc.ICs4events.gait(idx_detection).startSearch = startSearch;
EEG.etc.ICs4events.gait(idx_detection).endSearch = endSearch;
EEG.etc.ICs4events.gait(idx_detection).medianFilterLength = medianFilterLength;
EEG.etc.ICs4events.gait(idx_detection).minStepDistance = minStepDistance;
EEG.etc.ICs4events.gait(idx_detection).minPeakDuration = minPeakDuration;
EEG.etc.ICs4events.gait(idx_detection).maxPeakDuration = maxPeakDuration;
EEG.etc.ICs4events.gait(idx_detection).minPeakProminence = minPeakProminence;

%% detect steps

% median smoothing to preserve sharp peaks
smoothICwalking = smoothdata(EEG.icaact(idx_gait,startSearch:endSearch),'movmedian',medianFilterLength*EEG.srate);

% search peaks with a minimum peak distance shortly below step length
[~, locs, w, p] = findpeaks(smoothICwalking,'MinPeakDistance',minStepDistance*EEG.srate,...
    'MinPeakWidth',minPeakDuration*EEG.srate,'MaxPeakWidth',maxPeakDuration*EEG.srate,'MinPeakProminence',minPeakProminence);
[~, locs_flipped, w_flipped, p_flipped] = findpeaks(-smoothICwalking,'MinPeakDistance',minStepDistance*EEG.srate,...
    'MinPeakWidth',minPeakDuration*EEG.srate,'MaxPeakWidth',maxPeakDuration*EEG.srate,'MinPeakProminence',minPeakProminence);

plothandles = figure('position',[111 365 1083 419],'color','w');
subplot(2,2,1:2);
if mean(w_flipped) < mean(w)
    
    % smaller peak width is the sharp peak we want to detect as onset
    locs = locs_flipped;
    w = w_flipped;
    p = p_flipped;
    % find flipped peaks again for plot
    findpeaks(-smoothICwalking,EEG.times(startSearch:endSearch)/1000,...
        'MinPeakDistance',minStepDistance,'MinPeakProminence',minPeakProminence,...
    'MinPeakWidth',minPeakDuration,'MaxPeakWidth',maxPeakDuration,'annotate','extents')
else
    % find original peaks again for plot
    findpeaks(smoothICwalking,EEG.times(startSearch:endSearch)/1000,...
        'MinPeakDistance',minStepDistance,'MinPeakProminence',minPeakProminence,...
    'MinPeakWidth',minPeakDuration,'MaxPeakWidth',maxPeakDuration,'annotate','extents')
end

title(['steps, IC# ' num2str(idx_gait)])
xlabel('time [s]')
ylabel('smoothed IC activity [\muV]')

% plot prominences and widths
subplot(223);
histogram(p,'binwidth',0.1)
title('all prominences')
hold on
l=plot([minPeakProminence minPeakProminence],ylim,'r');
legend(l,'min prominence')
xlabel('\muV')

subplot(224);
histogram(w,'binwidth',2)
title('all widths')
hold on
l(1)=plot([minPeakDuration minPeakDuration],ylim,'r');
l(2)=plot([maxPeakDuration maxPeakDuration],ylim,'k');
legend(l,{'min width','max width'})
xlabel('seconds')

% store
EEG.etc.ICs4events.gait(idx_detection).prominences = p;
EEG.etc.ICs4events.gait(idx_detection).widths = w;

% add offset of start search area
locs = locs + (startSearch-1);

%% create events for the gait-related activity

EEG.etc.ICs4events.gait(idx_detection).latencies = locs;

for latency = locs
    i = numel(EEG.event) + 1;
    EEG.event(i).type = ['step:' num2str(idx_gait) ];
    EEG.event(i).latency = latency;
    EEG.event(i).duration = 1/EEG.srate;
end      

EEG = eeg_checkset(EEG, 'eventconsistency');

%% plot final steps only

rejcomps = 1:size(EEG.icaact,1);
rejcomps([idx_gait]) = [];

EEG_plot = pop_subcomp( EEG, rejcomps, 0);
pop_eegplot( EEG_plot, 0, 1, 1); 
plothandles(end+1) = gcf;

disp('Step detection done!')