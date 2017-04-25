% bemobil_preprocess() - Preprocessing of EEG data
%
% Usage:
%   >>  EEG = bemobil_preprocess( EEG, 'key', 'value');
%
% Inputs:
%   channel_locations_file     - channel_locations file (with path)
%   channels_to_remove    - e.g. with the current (22/11/2016) BeMoBIL layout {'N29' 'N30' 'N31'} are not in use
%   eog_channels          - ['G16' 'G32']
%   locutoff:             - low cut off frequency for firfilt filering, if NaN, no filter will be applied
%   highcutoff            - high cut of frequency, if NaN, no filter will be applied
%   resample_freq         - Resample frequency (Hz), if NaN or nonexistent, no resampling
%   out_filename          - output filename, if none, 'preprocessed.set' will be used
%   out_filepath          - output filepath for the merged file (OPTIONAL ARGUMENT - DEFAULT VALUE IS SAME FOLDER)
%
% Outputs:
%   EEG     - EEGLAB EEG structure
%
% See also:
%   POP_BEMOBIL_PREPROCESS, EEGLAB
% 
% Authors: Lukas Gehrke, Marius Klug, Friederike Hohlefeld, 2017

function [ ALLEEG EEG CURRENTSET ] = bemobil_preprocess(ALLEEG, EEG, CURRENTSET, channel_locations_filepath, channels_to_remove, eog_channels, locutoff, highcutoff, resample_freq, out_filename, out_filepath)

if nargin < 1
    help bemobil_preprocess;
    return;
end;

if ~exist('out_filename', 'var') out_filename = 'preprocessed.set'; end;
if ~exist('out_filepath', 'var') out_filepath = EEG.filepath; end;

% make sure output folder exists, nothing changes, if yes
mkdir(out_filepath);

% check if preprocessed file already exist and break if it does

dir_files = dir(out_filepath);
if ismember(out_filename, {dir_files.name})
    error(['Warning: ' out_filename ' file already exists in: ' out_filepath '. ' 'Exiting...']);
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
if ~isempty(channel_locations_filepath)
    EEG = pop_chanedit(EEG, 'load',...
        {channel_locations_filepath 'filetype' 'autodetect'});
    disp('Imported channel locations.');
    EEG.urchanlocs = EEG.chanlocs;
end

% 1d) change channel names in standard MoBI montage declaring the EOG
% channels
if ~isempty(eog_channels)
    for n = 1:length(EEG.chanlocs)
        if ismember(lower(EEG.chanlocs(n).labels), lower(eog_channels))
            EEG.chanlocs(n).type = strcat('EOG');
            disp(['Changed channel type: ', EEG.chanlocs(n).labels, ' to EOG electrode.']);
        end
    end
    EEG = eeg_checkset( EEG );
end

% 2. Resample/downsample to 250Hz if no ther resampling frequency is
% provided
if exist('resample_freq','var') && ~isempty(resample_freq)
    EEG = pop_resample(EEG, resample_freq);
    EEG = eeg_checkset( EEG );
    disp(['Resampled data to: ', num2str(resample_freq), 'Hz.']);
end

% 3. Filtering
% Resources https://sccn.ucsd.edu/wiki/Firfilt_FAQ

% highpass

if ~isempty(locutoff)
   
    figure;
    [EEG, com, b] = pop_eegfiltnew(EEG, locutoff, 0, [], 0, [], 1);
    EEG = eeg_checkset( EEG );
    saveas(gcf,[out_filepath '\preprocessing_filter_response_highpass_' num2str(locutoff) 'Hz']);
    
    split1 = strsplit(com, ' ');
    split2 = strsplit(split1{6}, ',');
    highpass_order = str2num(split2{1});
    highpass_cutoff = locutoff/2; % in eeglab the specified cutoff is the passband edge
    highpass_passband = locutoff;
    highpass_transition_bandwidth = locutoff;
    
    disp(['Highpass filtered the data with ' num2str(highpass_cutoff) 'Hz cutoff, '...
        num2str(highpass_transition_bandwidth) 'Hz transition bandwidth, '...
        num2str(highpass_passband) 'Hz passband edge, and '...
        num2str(highpass_order) ' order.']);
    
    EEG.etc.filter.highpass.cutoff = highpass_cutoff;
    EEG.etc.filter.highpass.transition_bandwidth = highpass_transition_bandwidth;
    EEG.etc.filter.highpass.passband = highpass_passband;
    EEG.etc.filter.highpass.order = highpass_order;
    close;
else
    
    EEG.etc.filter.highpass = 'not applied';
    
end


% lowpass

if ~isempty(highcutoff)
    
    if highcutoff > (EEG.srate/2) - 1
    disp('Warning: Cannot filter higher than Nyquist frequency.');
    highcutoff = (EEG.srate/2) - 1;
    disp(['Now continuing with highest possible frequency: ' num2str(highcutoff)]);
    end
   
    figure;
    [EEG, com, b] = pop_eegfiltnew(EEG, 0, highcutoff, [], 0, [], 1);
    EEG = eeg_checkset( EEG );
    saveas(gcf,[out_filepath '\preprocessing_filter_response_lowpass_' num2str(highcutoff) 'Hz']);
    
    split1 = strsplit(com, ' ');
    split2 = strsplit(split1{6}, ',');
    lowpass_order = str2num(split2{1});
    lowpass_transition_bandwidth = highcutoff*0.25/2;
    lowpass_cutoff = highcutoff + lowpass_transition_bandwidth; % in eeglab the specified cutoff is the passband edge
    lowpass_passband = highcutoff;
    
    disp(['Highpass filtered the data with ' num2str(lowpass_cutoff) 'Hz cutoff, '...
        num2str(lowpass_transition_bandwidth) 'Hz transition bandwidth, '...
        num2str(lowpass_passband) 'Hz passband edge, and '...
        num2str(lowpass_order) ' order.']);
    
    EEG.etc.filter.lowpass.cutoff = lowpass_cutoff;
    EEG.etc.filter.lowpass.transition_bandwidth = lowpass_transition_bandwidth;
    EEG.etc.filter.lowpass.passband = lowpass_passband;
    EEG.etc.filter.lowpass.order = lowpass_order;
    close;
else
    
    EEG.etc.filter.lowpass = 'not applied';
    
end

EEG.etc.filter.type = 'Hamming windowed sinc FIR filter';

%save data and stop function so manual channel rejection is possible
[ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, CURRENTSET, 'gui', 'off');
EEG = eeg_checkset( EEG );
EEG = pop_saveset( EEG, 'filename',out_filename,'filepath', out_filepath);
disp('...done');
[ALLEEG EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);
end

