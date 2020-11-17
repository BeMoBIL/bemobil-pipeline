% bemobil_filter() - Filtering of EEG data. Information about the filter is
% stired in EEG.etc.filter. Filter order can optionally be specified,
% otherwise EEGLAB default filters are being used. The filters are being
% applied successively to not use a band-pass filter with the same order.
% First, the highpass, then the lowpass filter are being run. A zero-phase
% Hamming window FIR filter is being used.
%
% Usage:
%   >>  [ ALLEEG EEG CURRENTSET ] = bemobil_filter(ALLEEG, EEG, CURRENTSET, lowerPassbandEdge, higherPassbandEdge);
%   >>  [ ALLEEG EEG CURRENTSET ] = bemobil_filter(ALLEEG, EEG, CURRENTSET, lowerPassbandEdge, higherPassbandEdge, out_filename, out_filepath);
%   >>  [ ALLEEG EEG CURRENTSET ] = bemobil_filter(ALLEEG, EEG, CURRENTSET, lowerPassbandEdge, higherPassbandEdge, out_filename, out_filepath, highPassFilterOrder, lowPassFilterOrder);
%
% Inputs:
%   ALLEEG                  - complete EEGLAB data set structure
%   EEG                     - current EEGLAB EEG structure
%   CURRENTSET              - index of current EEGLAB EEG structure within ALLEEG
%   lowerPassbandEdge       - low cut off frequency for firfilt filering, if [], no filter will be applied
%   higherPassbandEdge      - high cut of frequency, if [], no filter will be applied
%   out_filename            - output filename (OPTIONAL ARGUMENT)
%   out_filepath            - output filepath (OPTIONAL ARGUMENT - File will only be saved on disk
%       if both a name and a path are provided and not empty)
%	highPassFilterOrder		- OPTIONAL specify the order of the highpassfilter
%								(e.g. 1650 is a transition bandwidth of 0.5Hz)
%	lowPassFilterOrder		- OPTIONAL specify the order of the lowpassfilter
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
% Authors: Marius Klug, 2020

function [ ALLEEG EEG CURRENTSET ] = bemobil_filter(ALLEEG, EEG, CURRENTSET, lowerPassbandEdge, higherPassbandEdge,...
    out_filename, out_filepath, highPassFilterOrder, lowPassFilterOrder)

% only save a file on disk if both a name and a path are provided
save_file_on_disk = (exist('out_filename', 'var') && exist('out_filepath', 'var')) && ~isempty(out_filename) && ~isempty(out_filepath);

% check if file already exist and show warning if it does
if save_file_on_disk
    mkdir(out_filepath); % make sure that folder exists, nothing happens if so
    dir_files = dir(out_filepath);
    if ismember(out_filename, {dir_files.name})
        warning([out_filename ' file already exists in: ' out_filepath '. File will be overwritten...']);
    end
end

% Resources https://sccn.ucsd.edu/wiki/Firfilt_FAQ

if ~isempty(lowerPassbandEdge) ||  ~isempty(higherPassbandEdge)
    EEG.etc.filter.type = 'Hamming windowed zero-phase FIR filter';
else
    error('No filter cutoffs specified, what was your plan here?!')
end

if ~exist('highPassFilterOrder','var')
	highPassFilterOrder = [];
end
if ~exist('lowPassFilterOrder','var')
	lowPassFilterOrder = [];
end

% highpass

if ~isempty(lowerPassbandEdge)
   
    figure;
    [EEG, com, b] = pop_eegfiltnew(EEG, lowerPassbandEdge, 0, highPassFilterOrder, 0, [], 1);
    EEG = eeg_checkset( EEG );
    
    if save_file_on_disk; saveas(gcf,[out_filepath '\filter_response_highpass']); end
    
    split1 = strsplit(com, ' ');
    split2 = strsplit(split1{4}, ',');
	highpass_order = str2num(split2{6});
    highpass_passband = lowerPassbandEdge;
	
	if isempty(highPassFilterOrder)
		maxDf = lowerPassbandEdge;
		highpass_transition_bandwidth = min([max([maxDf * 0.25 2]) maxDf]); % this comes from pop_eegfiltnew
	else
		highpass_transition_bandwidth = 3.3 / highpass_order * EEG.srate; % this comes from pop_eegfiltnew
	end
	
    highpass_cutoff = highpass_passband-highpass_transition_bandwidth/2;
    
    disp(['Highpass filtered the data with ' num2str(highpass_passband) 'Hz passband edge, '...
        num2str(highpass_transition_bandwidth) 'Hz transition bandwidth, '...
        num2str(highpass_cutoff) 'Hz cutoff, and '...
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
        % removing and remaking the field is necessary for the order of the struct fields to be identical
		% yes I know, but I just have to.
        temp = EEG.etc.filter.highpass;
        EEG.etc.filter = rmfield(EEG.etc.filter, 'highpass');
        EEG.etc.filter.highpass = temp;
    end
    
end


% lowpass

if ~isempty(higherPassbandEdge)
    
    if higherPassbandEdge > (EEG.srate/2) - 1
    disp('Warning: Cannot filter higher than Nyquist frequency.');
    higherPassbandEdge = (EEG.srate/2) - 1;
    disp(['Now continuing with highest possible frequency: ' num2str(higherPassbandEdge)]);
    end
   
    figure;
    [EEG, com, b] = pop_eegfiltnew(EEG, 0, higherPassbandEdge, lowPassFilterOrder, 0, [], 1);
    EEG = eeg_checkset( EEG );
    if save_file_on_disk; saveas(gcf,[out_filepath '\filter_response_lowpass']); end
	
	split1 = strsplit(com, ' ');
    split2 = strsplit(split1{4}, ',');
	lowpass_order = str2num(split2{6});
    lowpass_passband = higherPassbandEdge;
	
	if isempty(lowPassFilterOrder)
		maxDf = EEG.srate/2 - higherPassbandEdge;
		lowpass_transition_bandwidth = min([max([higherPassbandEdge * 0.25 2]) maxDf]); % this comes from pop_eegfiltnew
	else
		lowpass_transition_bandwidth = 3.3 / lowpass_order * EEG.srate; % this comes from pop_eegfiltnew
	end
	
    lowpass_cutoff = higherPassbandEdge + lowpass_transition_bandwidth/2; % in eeglab the specified cutoff is the passband edge
    
    disp(['Lowpass filtered the data with ' num2str(lowpass_passband) 'Hz passband edge, ' ...
        num2str(lowpass_transition_bandwidth) 'Hz transition bandwidth, '...
        num2str(lowpass_cutoff) 'Hz cutoff, and '...
        num2str(lowpass_order) ' order.']);
    
    % removing and remaking the filter struct field is necessary for the order of the struct fields to be identical
	% yes I know, but I just have to.
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
