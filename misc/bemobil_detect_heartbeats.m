% bemobil_detect_heartbeats detects heart beats from an ECG channel using Pan-Tompkins algorithm and adds them as events
% to the data set.
%
% The implementation is a modified version of an implementation avaliable at
% https://www.researchgate.net/publication/313673153_Matlab_Implementation_of_pan_tompkins_ECG_QRS_detector
% References :
% [1] Sedghamiz. H, "Matlab Implementation of Pan Tompkins ECG QRS detector.",2014. (See researchgate)
% [2] PAN.J, TOMPKINS. W.J,"A Real-Time QRS Detection Algorithm" IEEE TRANSACTIONS ON BIOMEDICAL ENGINEERING, VOL. BME-32, NO. 3, MARCH 1985.
%
% Input:
%       EEG_heart           - EEGLAB dataset containing ECG data
%       idx_heart           - OPTIONAL channel index of ECG (default = 1)
%       highpass            - OPTIONAL highpass filter (default = 1)
%       lowpass             - OPTIONAL lowpass filter (default = 40)
%       createplots         - OPTIONAL boolean whether plots should be created (default = 1)
%
% Output:
%       EEG_heart           - EEGLAB dataset containing detected "heartbeat" events
%       plothandles         - handle to the three created plots to allow script-based saving and closing 
%
% SEE ALSO:
%   pan_tompkins


function [EEG_heart, plothandles] = bemobil_detect_heartbeats(EEG_heart,idx_heart,highpass,lowpass,createplots)

if nargin == 0
    help bemobil_ECG_analysis
    return
end

if ~isfield(EEG_heart,'data')
    error('Wrong data type for EEG_heart input! Needs to be an EEGLAB struct.')
end

if ~exist('idx_heart','var') || ~isvector(idx_heart)
    disp('Assuming only one channel exists in the data.')
    idx_heart = 1;
end

if ~exist('highpass','var') || ~isvector(highpass)
    highpass = 1;
    disp(['Using a ' num2str(highpass) ' Hz highpass filter.'])
end

if ~exist('lowpass','var') || ~isvector(lowpass)
    lowpass = 40;
    disp(['Using a ' num2str(lowpass) ' Hz lowpass filter.'])
end

if ~exist('createplots','var') || ~isscalar(createplots)
    createplots = 1;
end

%% use Pan-Tompkins algorithm

disp(['Detecting heart beats of channel: [' num2str(idx_heart) '], with label: ' EEG_heart.chanlocs(idx_heart).labels...
    ' using the Pan-Tompkins algorithm'])
disp(' ')
disp('--------------------------- PLEASE CITE ----------------------------')
disp(' ')
disp('[1] Sedghamiz. H, "Matlab Implementation of Pan Tompkins ECG QRS detector.",2014. (See researchgate)')
disp('[2] PAN.J, TOMPKINS. W.J,"A Real-Time QRS Detection Algorithm" IEEE TRANSACTIONS ON BIOMEDICAL ENGINEERING, VOL. BME-32, NO. 3, MARCH 1985.')
disp(' ')
disp('--------------------------- PLEASE CITE ----------------------------')
disp(' ')

[~,qrs_i_raw,~,plothandles] = pan_tompkins(EEG_heart.data(idx_heart,:),EEG_heart.srate,createplots,highpass,lowpass);

%% add events in EEG_heart set

disp('Adding "heartbeat" events to data set.')

for i_heart = 1:length(qrs_i_raw)
    
    i = numel(EEG_heart.event) + 1;
    EEG_heart.event(i).type = 'heartbeat';
    EEG_heart.event(i).latency = qrs_i_raw(i_heart);
    EEG_heart.event(i).duration = 1/EEG_heart.srate;
    
end

EEG_heart = eeg_checkset(EEG_heart, 'eventconsistency');

% store info
EEG_heart.etc.heartbeat_detection.idx_heart = idx_heart;
EEG_heart.etc.heartbeat_detection.highpass = highpass;
EEG_heart.etc.heartbeat_detection.lowpass = lowpass;

disp('Heart beat detection done!')

