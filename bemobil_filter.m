% bemobil_filter() - Filtering of EEG data
%
% Usage:
%   >>  [ ALLEEG EEG CURRENTSET ] = bemobil_filter(ALLEEG, EEG, CURRENTSET, locutoff, highcutoff);
%   >>  [ ALLEEG EEG CURRENTSET ] = bemobil_filter(ALLEEG, EEG, CURRENTSET, locutoff, highcutoff, out_filename, out_filepath);
%
% Inputs:
%   ALLEEG                  - complete EEGLAB data set structure
%   EEG                     - current EEGLAB EEG structure
%   CURRENTSET              - index of current EEGLAB EEG structure within ALLEEG
%   lowcutoff               - low cut off frequency for firfilt filering, if [], no filter will be applied
%   highcutoff              - high cut of frequency, if [], no filter will be applied
%   out_filename            - output filename (OPTIONAL ARGUMENT)
%   out_filepath            - output filepath (OPTIONAL ARGUMENT - File will only be saved on disk
%       if both a name and a path are provided)
%
% Outputs:
%   ALLEEG                  - complete EEGLAB data set structure
%   EEG                     - current EEGLAB EEG structure
%   Currentset              - index of current EEGLAB EEG structure within ALLEEG
%
%   .set data file of current EEGLAB EEG structure stored on disk (OPTIONALLY)
%
% See also:
%   EEGLAB, pop_eegfiltnew, bemobil_preprocess
% 
% Authors: Marius Klug, 2017

function [ ALLEEG EEG CURRENTSET ] = bemobil_filter(ALLEEG, EEG, CURRENTSET, lowcutoff, highcutoff, out_filename, out_filepath)

% only save a file on disk if both a name and a path are provided
save_file_on_disk = (exist('out_filename', 'var') && exist('out_filepath', 'var'));

% check if file already exist and show warning if it does
if save_file_on_disk
    mkdir(out_filepath); % make sure that folder exists, nothing happens if so
    dir_files = dir(out_filepath);
    if ismember(out_filename, {dir_files.name})
        warning([out_filename ' file already exists in: ' out_filepath '. File will be overwritten...']);
    end
end

% Resources https://sccn.ucsd.edu/wiki/Firfilt_FAQ

% highpass

if ~isempty(lowcutoff)
   
    figure;
    [EEG, com, b] = pop_eegfiltnew(EEG, lowcutoff, 0, [], 0, [], 1);
    EEG = eeg_checkset( EEG );
    saveas(gcf,[out_filepath '\preprocessing_filter_response_highpass_' num2str(lowcutoff) 'Hz']);
    
    split1 = strsplit(com, ' ');
    split2 = strsplit(split1{6}, ',');
    highpass_order = str2num(split2{1});
    highpass_cutoff = lowcutoff/2; % in eeglab the specified cutoff is the passband edge
    highpass_passband = lowcutoff;
    highpass_transition_bandwidth = lowcutoff;
    
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

% new data set in EEGLAB
[ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, CURRENTSET, 'gui', 'off');
EEG = eeg_checkset( EEG );

% save on disk
if save_file_on_disk
    EEG = pop_saveset( EEG, 'filename',out_filename,'filepath', out_filepath);
    disp('...done');
end

[ALLEEG EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);
