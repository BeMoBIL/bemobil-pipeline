% bemobil_unmix_electrode_mixups() - Rearranges channels that have
% been mixed up physically, during location digitization, or as data
% signal. If only location was mixed up during digitization, leave data
% matrix empty, if only data cables have been mixed up, leave location
% empty, if electrodes have been physically misplaced on the scalp but have
% correct location, use both data and location matrices.
%
% Usage:
%   >>  [ALLEEG, EEG, CURRENTSET] = bemobil_unmix_electrode_mixups( EEG , ALLEEG, CURRENTSET, data_unmix_matrix, location_unmix_matrix)
%   >>  [ALLEEG, EEG, CURRENTSET] = bemobil_unmix_electrode_mixups( EEG , ALLEEG, CURRENTSET, data_unmix_matrix, location_unmix_matrix, out_filename, out_filepath)
%
% Inputs:
%   ALLEEG                  - complete EEGLAB data set structure
%   EEG                     - current EEGLAB EEG structure
%   CURRENTSET              - index of current EEGLAB EEG structure within ALLEEG
%   data_unmix_matrix       - matrix of channels which should de-mix DATA
%       (e.g. [12 16 24; 16 24 12] meaning that  channel 12 will get the data
%       of channel 16, channel 16 the data of channel 24, and channel 24 the 
%       data of channel 12)
%   location_unmix_matrix   - matrix of channels which should de-mix
%       LOCATION (e.g. [12 16 24; 16 24 12] meaning that  channel 12 will get
%       the location of channel 16, channel 16 the location of channel 24, 
%       and channel 24 the location of channel 12)
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
% Authors: Marius Klug, 2018

function [ALLEEG, EEG, CURRENTSET] = bemobil_unmix_electrode_mixups( EEG , ALLEEG, CURRENTSET, data_unmix_matrix, location_unmix_matrix, out_filename, out_filepath)

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

if ~isempty(EEG.icaweights)

    error('ICA weights present. First de-mix channels and then calculate ICA!');
    
end

newEEG = EEG;
    
if ~isempty(data_unmix_matrix)
    % this is necessary if data cables have been wrongly plugged in
    
    newEEG.data(data_unmix_matrix(1,:),:) = EEG.data(data_unmix_matrix(2,:),:);
    
    newEEG.etc.demixed_channels.data = data_unmix_matrix;
    
end

if ~isempty(location_unmix_matrix)
    % this is necessary if the location digitized by the vicra xensor
    % system has been mixed up
    
    for channel = 1:size(location_unmix_matrix,2)
       
        newEEG.chanlocs(location_unmix_matrix(1,channel)).X             = EEG.chanlocs(location_unmix_matrix(2,channel)).X;
        newEEG.chanlocs(location_unmix_matrix(1,channel)).Y             = EEG.chanlocs(location_unmix_matrix(2,channel)).Y;
        newEEG.chanlocs(location_unmix_matrix(1,channel)).Z             = EEG.chanlocs(location_unmix_matrix(2,channel)).Z;
        newEEG.chanlocs(location_unmix_matrix(1,channel)).sph_phi       = EEG.chanlocs(location_unmix_matrix(2,channel)).sph_phi;
        newEEG.chanlocs(location_unmix_matrix(1,channel)).sph_radius    = EEG.chanlocs(location_unmix_matrix(2,channel)).sph_radius;
        newEEG.chanlocs(location_unmix_matrix(1,channel)).theta         = EEG.chanlocs(location_unmix_matrix(2,channel)).theta;
        newEEG.chanlocs(location_unmix_matrix(1,channel)).radius        = EEG.chanlocs(location_unmix_matrix(2,channel)).radius;
        newEEG.chanlocs(location_unmix_matrix(1,channel)).sph_theta     = EEG.chanlocs(location_unmix_matrix(2,channel)).sph_theta;
        
    end
    
    newEEG.etc.demixed_channels.locations = location_unmix_matrix;
end

% new data set in EEGLAB
[ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, newEEG, CURRENTSET, 'gui', 'off');
EEG = eeg_checkset( EEG );

% save on disk
if save_file_on_disk
    EEG = pop_saveset( EEG, 'filename',out_filename,'filepath', out_filepath);
    disp('...done');
end

[ALLEEG EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);