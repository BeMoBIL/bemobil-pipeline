% bemobil_interp_avref_copy_spatial_filter() - Interpolates missing channels with spherical interpolation,
% average references the data and, if given, copies the ICA weights from a given data set into the
% current EEG data set (either both data sets have to have no de-mixed channels before or both
% have to have it. Requirements (otherwise the ICA weights will either not work or work but
% be incorrect): both data sets have to be at the same stage of processing, meaning:
%
% - Same reference
% - Same chanloc interpolations
% - Same channel data interpolations
% - Same de-mixed channels (in case something was mixed up before)
%
% The weights might have been created by another algorithm than ICA, but due to the EEGLAB structure
% they must be stored as if the were ICA weights. In EEG.etc.spatial_filter all the necessary
% information about the algorithm is stored (is also copied).
%
% Usage:
%   >>  [ALLEEG, EEG, CURRENTSET] = bemobil_interp_avref_copyICA( EEG , ALLEEG, CURRENTSET, channels_to_interpolate)
%   >>  [ALLEEG, EEG, CURRENTSET] = bemobil_interp_avref_copyICA( EEG , ALLEEG, CURRENTSET, channels_to_interpolate, EEG_set_to_copy_spatial_filter, copy_dipfit, copy_reject)
%   >>  [ALLEEG, EEG, CURRENTSET] = bemobil_interp_avref_copyICA( EEG , ALLEEG, CURRENTSET, channels_to_interpolate, EEG_set_to_copy_spatial_filter, copy_dipfit, copy_reject, out_filename, out_filepath)
%
% Inputs:
%   ALLEEG                  - complete EEGLAB data set structure
%   EEG                     - current EEGLAB EEG structure
%   CURRENTSET              - index of current EEGLAB EEG structure within ALLEEG
%   EEG_set_to_copy_ICA     - data set that should donate the ICA weights; OR [], then no ICA weights will be copied
%   channels_to_interpolate - vector of channel numbers that should be interpolated; if [],
%       attempts to interpolate all missing (already deleted) channels from urchanlocs
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
%   EEGLAB, bemobil_copy_spatial_filter, pop_interp, pop_reref, pop_interp, 
% 
% Authors: Lukas Gehrke, Marius Klug, 2020

function [ALLEEG, EEG, CURRENTSET] = bemobil_interp_avref_copy_spatial_filter( EEG , ALLEEG, CURRENTSET, channels_to_interpolate, EEG_set_to_copy_spatial_filter, copy_dipfit, copy_reject, out_filename, out_filepath)

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

% Interpolate channels with spherical interpolation, reference to average of EEG channels (ignore EOG)
[ALLEEG, EEG, CURRENTSET] = bemobil_interp_avref( EEG , ALLEEG, CURRENTSET, channels_to_interpolate);

% Copy spatial filter weights
if ~isempty(EEG_set_to_copy_spatial_filter) && ~isempty(copy_dipfit) && ~isempty(copy_reject)
    
    [ALLEEG, EEG, CURRENTSET] = bemobil_copy_spatial_filter( EEG , ALLEEG, CURRENTSET, EEG_set_to_copy_spatial_filter, copy_dipfit, copy_reject);
    

else
    disp('No data set to copy ICA weights from or no booleans if copy dipfit and rejections is provided. Skipping this step.');
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
