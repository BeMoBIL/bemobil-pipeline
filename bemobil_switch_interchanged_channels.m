% bemobil_switch_interchanged_channels() - Rearranges channels that have
% been switched in location or data signal
%
% Usage:
%   >>  [ALLEEG, EEG, CURRENTSET] = bemobil_switch_interchanged_channels( EEG , ALLEEG, CURRENTSET, channelsNew, channelsOriginal, process)
%   >>  [ALLEEG, EEG, CURRENTSET] = bemobil_switch_interchanged_channels( EEG , ALLEEG, CURRENTSET, channelsNew, channelsOriginal, process, out_filename, out_filepath)
%
% Inputs:
%   ALLEEG                  - complete EEGLAB data set structure
%   EEG                     - current EEGLAB EEG structure
%   CURRENTSET              - index of current EEGLAB EEG structure within ALLEEG
%   channelsNew             - vector of channel numbers how they should be correctly (e.g. [12 16])
%   chanelsOriginal         - vectors of the channels numbers how they are now (e.g. [16 12])
%       the original channels will be projected to the new corresponding channels (e.g. 16 -> 12 and
%       12 -> 16)
%   process                 - which process should be used for the switched
%       channels (can be either 'exchange labels, type & channel number' 
%       for physically misplaced electrodes, 'exchange data' for wrongly 
%       plugged in data cables, or 'exchange locations' for mix-up in the
%       electrode location digitization)
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
%   EEGLAB
% 
% Authors: Marius Klug, 2017

function [ALLEEG, EEG, CURRENTSET] = bemobil_switch_interchanged_channels( EEG , ALLEEG, CURRENTSET, channelsNew, channelsOriginal, process, out_filename, out_filepath)

error('Deprecated version, use ''bemobil_unmix_electrode_mixups()'' istead of ''bemobil_switch_interchanged_channels()''. If you know what you''re doing, uncomment this line.')

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

% check if all new channels have a corresponding old channel

if length(channelsNew)~=length(channelsOriginal)
    error('There must be the same number of new channels and their corresponding original channels!');
end

if ~isempty(EEG.icaweights)

    error('ICA weights present. First interchange channels and then calculate ICA!');
    
end

newEEG = EEG;
newEEG.etc.interchanged_channels.process = process;

if strcmp(process, 'exchange labels, type & channel number')
    % this is neccessary if e.g. electrodes have been physically misplaced, 
    % G16 (eye) is at the position of G10 (forehead) for example
    
    tempNewEEG = newEEG;
    

    % first interchange labels and type, then the data is connected to the
    % correct labes/type
    
    for channel = 1:length(channelsNew)
       
        tempNewEEG.chanlocs(channelsNew(channel)).labels = EEG.chanlocs(channelsOriginal(channel)).labels;
        tempNewEEG.chanlocs(channelsNew(channel)).type = EEG.chanlocs(channelsOriginal(channel)).type;
        
    end
    
    % then interchange data & the complete chanloc info, this accounts for
    % the fact that the channel number has to be interchanged back such
    % that everything is in the correct order (otherwise e.g. channel g16,
    % which was g10, still would be the 10th data channel)
    newEEG.chanlocs(channelsNew) = tempNewEEG.chanlocs(channelsOriginal);
    newEEG.data(channelsNew,:) = tempNewEEG.data(channelsOriginal,:);
    
elseif strcmp(process, 'exchange data')
    % this is necessary if data cables have been wrongly plugged in
    
    newEEG.data(channelsNew,:) = EEG.data(channelsOriginal,:);
    
elseif strcmp(process, 'exchange locations')
    % this is necessary if the location digitized by the vicra xensor
    % system has been mixed up
    
    for channel = 1:length(channelsNew)
       
        newEEG.chanlocs(channelsNew(channel)).X = EEG.chanlocs(channelsOriginal(channel)).X;
        newEEG.chanlocs(channelsNew(channel)).Y = EEG.chanlocs(channelsOriginal(channel)).Y;
        newEEG.chanlocs(channelsNew(channel)).Z = EEG.chanlocs(channelsOriginal(channel)).Z;
        newEEG.chanlocs(channelsNew(channel)).sph_phi = EEG.chanlocs(channelsOriginal(channel)).sph_phi;
        newEEG.chanlocs(channelsNew(channel)).sph_radius = EEG.chanlocs(channelsOriginal(channel)).sph_radius;
        newEEG.chanlocs(channelsNew(channel)).theta = EEG.chanlocs(channelsOriginal(channel)).theta;
        newEEG.chanlocs(channelsNew(channel)).radius = EEG.chanlocs(channelsOriginal(channel)).radius;
        newEEG.chanlocs(channelsNew(channel)).sph_theta = EEG.chanlocs(channelsOriginal(channel)).sph_theta;
        
    end
    
else
    
    error('U done fucked up. ''Process'' must be either ''exchange labels, type & channel number'', ''exchange data'', or ''exchange locations''!')
    
end

newEEG.etc.interchanged_channels.channels_original = channelsOriginal;
newEEG.etc.interchanged_channels.channels_new = channelsNew;

% new data set in EEGLAB
[ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, newEEG, CURRENTSET, 'gui', 'off');
EEG = eeg_checkset( EEG );

% save on disk
if save_file_on_disk
    EEG = pop_saveset( EEG, 'filename',out_filename,'filepath', out_filepath);
    disp('...done');
end

[ALLEEG EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);