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

if ~isempty(lowcutoff) ||  ~isempty(highcutoff)
    EEG.etc.filter.type = 'Hamming windowed sinc FIR filter (zero-phase)';
else
    error('No filter cutoffs specified, what was your plan here?!')
end

% highpass

if ~isempty(lowcutoff)
   
    figure;
    [EEG, com, b] = pop_eegfiltnew(EEG, lowcutoff, 0, [], 0, [], 1);
    EEG = eeg_checkset( EEG );
    
    if save_file_on_disk; saveas(gcf,[out_filepath '\filter_response_highpass']); end
    
    split1 = strsplit(com, ' ');
    split2 = strsplit(split1{4}, ',');
    highpass_order = str2num(split2{3}) + 1;
    highpass_cutoff = lowcutoff/2; % in eeglab the specified cutoff is the passband edge
    highpass_passband = lowcutoff;
    highpass_transition_bandwidth = lowcutoff;
    
    disp(['Highpass filtered the data with ' num2str(highpass_cutoff) 'Hz cutoff, '...
        num2str(highpass_transition_bandwidth) 'Hz transition bandwidth, '...
        num2str(highpass_passband) 'Hz passband edge, and '...
        num2str(highpass_order) ' order.']);
    
    % removing and remaking the filed is necessary for the order of the struct fields to be identical
    if isfield(EEG.etc.filter,'highpass');  EEG.etc.filter = rmfield(EEG.etc.filter, 'highpass'); end 
    EEG.etc.filter.highpass.cutoff = highpass_cutoff;
    EEG.etc.filter.highpass.transition_bandwidth = highpass_transition_bandwidth;
    EEG.etc.filter.highpass.passband = highpass_passband;
    EEG.etc.filter.highpass.order = highpass_order;
    close;
else
    
    if ~isfield(EEG.etc.filter,'highpass')
        EEG.etc.filter.highpass = 'not applied'; 
    else
        % removing and remaking the filed is necessary for the order of the struct fields to be identical
        temp = EEG.etc.filter.highpass;
        EEG.etc.filter = rmfield(EEG.etc.filter, 'highpass');
        EEG.etc.filter.highpass = temp;
    end
    
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
    if save_file_on_disk; saveas(gcf,[out_filepath '\filter_response_lowpass']); end
    
    split1 = strsplit(com, ' ');
    split2 = strsplit(split1{4}, ',');
    lowpass_order = str2num(split2{3}) + 1;
    lowpass_transition_bandwidth = highcutoff*0.25;
    lowpass_cutoff = highcutoff + lowpass_transition_bandwidth/2; % in eeglab the specified cutoff is the passband edge
    lowpass_passband = highcutoff;
    
    disp(['Lowpass filtered the data with ' num2str(lowpass_cutoff) 'Hz cutoff, '...
        num2str(lowpass_transition_bandwidth) 'Hz transition bandwidth, '...
        num2str(lowpass_passband) 'Hz passband edge, and '...
        num2str(lowpass_order) ' order.']);
    
    % removing and remaking the filter struct field is necessary for the order of the struct fields to be identical
    if isfield(EEG.etc.filter,'lowpass'); EEG.etc.filter = rmfield(EEG.etc.filter, 'lowpass'); end
    EEG.etc.filter.lowpass.cutoff = lowpass_cutoff;
    EEG.etc.filter.lowpass.transition_bandwidth = lowpass_transition_bandwidth;
    EEG.etc.filter.lowpass.passband = lowpass_passband;
    EEG.etc.filter.lowpass.order = lowpass_order;
    close;
else
    
    if ~isfield(EEG.etc.filter,'lowpass')
        EEG.etc.filter.lowpass = 'not applied'; 
    else
        % removing and remaking the filed is necessary for the order of the struct fields to be identical
        temp = EEG.etc.filter.lowpass;
        EEG.etc.filter = rmfield(EEG.etc.filter, 'lowpass');
        EEG.etc.filter.lowpass = temp;
    end
    
end

% new data set in EEGLAB
[ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, CURRENTSET, 'gui', 'off');
EEG = eeg_checkset( EEG );

% save on disk
if save_file_on_disk
    EEG = pop_saveset( EEG, 'filename',out_filename,'filepath', out_filepath);
    disp('...done');
end

[ALLEEG EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);
