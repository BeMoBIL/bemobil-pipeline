% bemobil_preprocess() - Preprocessing of EEG data
%
% Usage:
%   >>  EEG = bemobil_preprocess( EEG, 'key', 'value');
%
% Inputs:
%   channel_locations     - filename of channel_locations file
%   channels_to_remove    - e.g. with the current (22/11/2016) BeMoBIL layout {'N29' 'N30' 'N31'} are not in use
%   eog_channels          - ['G16' 'G32']
%   locutoff:             - low cut off frequency for firfilt filering. If none, locutoff is set to 1
%   highcutoff            - high cut of frequency, if none, data will only be high pass filtered
%   resample_freq         - Resample frequency (Hz)
%    
% Outputs:
%   EEG     - EEGLAB EEG structure
%
% See also: 
%   POP_BEMOBIL_PREPROCESS, EEGLAB

function [ EEG ] = bemobil_preprocess(EEG, channel_locations, channels_to_remove, eog_channels, locutoff, highcutoff, resample_freq)

if nargin < 1
	help bemobil_preprocess;
	return;
end;	

% check if preprocessed file already exist and break if it does
out_filename = ['Preprocessed_' EEG.filename];
dir_files = dir(EEG.filepath);
if ismember(out_filename, {dir_files.name})
    error(['Warning: preprocessed file already exists in: ' EEG.filepath '. ' 'Exiting...']);
    %return; use only if warning is provided only on console with disp
end

% 1a) fill/copy all ur_structures with raw data (e.g. copy event to urevent)
EEG = eeg_checkset(EEG, 'makeur');

% 1b) remove unused neck electrodes from file (if BeMoBIL layout is used as is)
if ismember(channels_to_remove, {EEG.chanlocs.labels})
    % b) remove not needed channels import "corrected" chanlocs file
    EEG = pop_select( EEG,'nochannel', channels_to_remove);
    EEG = eeg_checkset( EEG );  
    disp(['Removed electrodes: ' channels_to_remove ' from the dataset.']);
end

% 1c) import chanlocs and copy to urchanlocs
if ~isempty(channel_locations)
    chanlocfile_path = [EEG.filepath channel_locations];
    if ~isempty(channel_locations)
        EEG = pop_chanedit(EEG, 'load',...
            {chanlocfile_path 'filetype' 'autodetect'},...
            'eval','chans = pop_chancenter( chans, [],[]);');
        disp('Imported channel locations and recentered.');
    end
    EEG.urchanlocs = EEG.chanlocs;
end

% 1d) change channel names in standard MoBI montage declaring the EOG
% channels
if ~isempty(channel_locations)
    for n = 1:length(EEG.chanlocs)            
        if ismember(EEG.chanlocs(n).labels, eog_channels)
            EEG.chanlocs(n).type = strcat('EOG');
            disp(['Changed channel type: ', EEG.chanlocs(n).labels, ' to EOG electrode./n']);
        end            
    end
    EEG = eeg_checkset( EEG );
end

% 2. Resample/downsample to 250Hz if no ther resampling frequency is
% provided
if isempty(resample_freq)
    resample_freq = 250;
end
EEG = pop_resample(EEG, resample_freq);
EEG = eeg_checkset( EEG );
disp(['Resampled data to: ', num2str(resample_freq), 'Hz.']);

% 3. Filtering 
% Resources https://sccn.ucsd.edu/wiki/Firfilt_FAQ
if isnan(locutoff)
    locutoff = 1;
end

if isnan(highcutoff)
    highcutoff = (resample_freq/2) - 1; % highest possible freq of interest due to Nyquist
elseif highcutoff > (resample_freq/2) - 1
    disp('Warning: Cannot filter higher than 2*Nyquist+1 frequency.');
    highcutoff = (resample_freq/2) - 1;
    disp(['Now continueing with highest possible frequency: ' num2str(highcutoff)]);
end

% plot filter response and save it
figure;
EEG = pop_eegfiltnew(EEG, locutoff, highcutoff, [], 0, [], 1);
EEG = eeg_checkset( EEG );
disp(['Filtered the data from ', num2str(locutoff), ' to ', num2str(highcutoff)]);
current_pwd = pwd;
cd(EEG.filepath);
saveas(gcf,[out_filename '_filter_response_' 'eegfiltnew_' num2str(locutoff) '-' num2str(highcutoff) '.eps'], 'psc2');
close;
cd(current_pwd);

%save data and stop function so manual channel rejection is possible
pop_saveset(EEG, strcat(EEG.filepath, out_filename));
end

