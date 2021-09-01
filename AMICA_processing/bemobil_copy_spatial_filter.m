% bemobil_copy_spatial_filter() - Copies the ICA weights, dipole fitting, and IC Label data from a given data set into
% another EEG data set. Requirements (otherwise the ICA weights will either not work or work but
% be incorrect): both data sets have to be at the same stage of processing, meaning:
%
% - Same reference
% - Same chanloc interpolations
% - Same channel data interpolations
% - Same interchanged channels (in case something was mixed up before)
%
% The weights might have been created by another algorithm than ICA, but due to the EEGLAB structure
% they must be stored as if the were ICA weights. In EEG.etc.spatial_filter all the necessary
% information about the algorithm is stored (is also copied).
%
% Usage:
%   >>  [ALLEEG, EEG, CURRENTSET] = bemobil_copy_spatial_filter( EEG , ALLEEG, CURRENTSET, EEG_set_to_copy_spatial_filter)
%   >>  [ALLEEG, EEG, CURRENTSET] = bemobil_copy_spatial_filter( EEG , ALLEEG, CURRENTSET, EEG_set_to_copy_spatial_filter, out_filename, out_filepath)
%
% Inputs:
%   ALLEEG                  - complete EEGLAB data set structure
%   EEG                     - current EEGLAB EEG structure
%   CURRENTSET              - index of current EEGLAB EEG structure within ALLEEG
%   EEG_set_to_copy_ICA     - data set that should donate the ICA weights
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
%   EEGLAB, bemobil_interp_avref_copyICA, pop_interp, pop_reref, pop_interp,
%
% Authors: Marius Klug, 2017

function [ALLEEG, EEG, CURRENTSET] = bemobil_copy_spatial_filter( EEG , ALLEEG, CURRENTSET, EEG_set_to_copy_spatial_filter, out_filename, out_filepath)

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

if ~isempty(EEG_set_to_copy_spatial_filter)
   
   % check if the requirements are met
   
   same_reference = isequal(EEG.etc.bemobil_reref,EEG_set_to_copy_spatial_filter.etc.bemobil_reref);
   same_interpolated_locations = true; % initialize as true, because they both might not have the fields, which is fine
   same_interpolated_channels = true;
   same_demixed_channels = true;
   
   % check the fields
   if isfield(EEG_set_to_copy_spatial_filter.etc,'interpolated_locations') % one has it
      if ~(isfield(EEG.etc,'interpolated_locations') && isequal(EEG.etc.interpolated_locations, EEG_set_to_copy_spatial_filter.etc.interpolated_locations)) % either the other doesn't have it or it's not the same
         same_interpolated_locations = false;
      end
   elseif isfield(EEG.etc,'interpolated_locations')
      if ~(isfield(EEG_set_to_copy_spatial_filter.etc,'interpolated_locations') && isequal(EEG.etc.interpolated_locations, EEG_set_to_copy_spatial_filter.etc.interpolated_locations))
         same_interpolated_locations = false;
      end
   end
   
   if isfield(EEG_set_to_copy_spatial_filter.etc,'interpolated_channels')
      if ~(isfield(EEG.etc,'interpolated_channels') && isequal(EEG.etc.interpolated_channels, EEG_set_to_copy_spatial_filter.etc.interpolated_channels))
         same_interpolated_channels = false;
      end
   elseif isfield(EEG.etc,'interpolated_channels')
      if ~(isfield(EEG_set_to_copy_spatial_filter.etc,'interpolated_channels') && isequal(EEG.etc.interpolated_channels, EEG_set_to_copy_spatial_filter.etc.interpolated_channels))
         same_interpolated_channels = false;
      end
   end
   
   if isfield(EEG_set_to_copy_spatial_filter.etc,'demixed_channels')
      if ~(isfield(EEG.etc,'demixed_channels') && isequal(EEG.etc.demixed_channels, EEG_set_to_copy_spatial_filter.etc.demixed_channels))
         same_demixed_channels = false;
      end
   elseif isfield(EEG.etc,'demixed_channels')
      if ~(isfield(EEG_set_to_copy_spatial_filter.etc,'demixed_channels') && isequal(EEG.etc.demixed_channels, EEG_set_to_copy_spatial_filter.etc.demixed_channels))
         same_demixed_channels = false;
      end
      
   end
   
   
   if ~all([same_reference,same_interpolated_locations,same_interpolated_channels,same_demixed_channels])
      
      error('Requirements not met. Type ''help bemobil_copy_spatial_filter'' for information!')
   end
   
   % requirements are met!
   
   disp('Copying ICA weights from provided data set.');
   EEG.icaweights = EEG_set_to_copy_spatial_filter.icaweights;
   EEG.icasphere = EEG_set_to_copy_spatial_filter.icasphere;
   
   % recompute the rest of ICA stuff
   EEG = eeg_checkset( EEG );
   
   if isfield(EEG_set_to_copy_spatial_filter.etc,'spatial_filter')
      
      EEG.etc.spatial_filter = EEG_set_to_copy_spatial_filter.etc.spatial_filter;
      
   end
   
   % copy dipfit stuff as well (if existing)
   
   if isfield(EEG_set_to_copy_spatial_filter,'dipfit')
      disp('Copying dipfit information into new dataset...')
      EEG.dipfit = EEG_set_to_copy_spatial_filter.dipfit;
      if isfield(EEG_set_to_copy_spatial_filter.etc,'dipfit')
         disp('Copying additional dipfit information (etc.dipfit) into new dataset...')
         EEG.etc.dipfit = EEG_set_to_copy_spatial_filter.etc.dipfit;
      else
         disp('Attempted to copy additional dipfit information in EEG.etc but they were not present!')
      end
   else
      disp('Attempted to copy dipfit information but they were not present.')
   end
   
   % copy component rejection information (if existing)
   if ~isempty(EEG_set_to_copy_spatial_filter.reject.gcompreject)
      disp('Copying IC rejection information into new dataset...')
      EEG.reject.gcompreject=EEG_set_to_copy_spatial_filter.reject.gcompreject;
   else
      disp('Attempted to copy component rejection information, but it was empty.')
   end
   
   % copy SASICA results (if existing)
   if isfield(EEG_set_to_copy_spatial_filter.reject,'SASICA')
      disp('Copying SASICA information into new dataset...')
      EEG.reject.SASICA=EEG_set_to_copy_spatial_filter.reject.SASICA;
   else
      disp('Attempted to copy SASICA information but they were not present.')
   end
   
   % copy ICLabel results (if existing)
   if isfield(EEG_set_to_copy_spatial_filter.etc,'ic_classification')
      disp('Copying ICLabel information into new dataset...')
      EEG.etc.ic_classification=EEG_set_to_copy_spatial_filter.etc.ic_classification;
   else
      disp('Attempted to copy ICLabel information but they were not present.')
   end
   
   % copy bad sample AMICA autoreject results (if existing)
   if isfield(EEG_set_to_copy_spatial_filter.etc,'bad_samples')
      disp('Copying bad samples information into new dataset...')
      EEG.etc.bad_samples=EEG_set_to_copy_spatial_filter.etc.bad_samples;
      EEG.etc.bad_samples_percent=EEG_set_to_copy_spatial_filter.etc.bad_samples_percent;
      EEG.etc.remove_data_intervals=EEG_set_to_copy_spatial_filter.etc.remove_data_intervals;
   else
      disp('Attempted to copy AMICA autoreject bad sample information but they were not present.')
   end
   
else
   error('No data set to copy ICA weights from is provided. Wat?!');
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

