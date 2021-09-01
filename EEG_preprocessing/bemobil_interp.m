% bemobil_interp() - Interpolates missing channels with spherical interpolation and
% average references the data.
%
% Usage:
%   >>  [ALLEEG, EEG, CURRENTSET] = bemobil_interp_avref( EEG , ALLEEG, CURRENTSET, channels_to_interpolate)
%   >>  [ALLEEG, EEG, CURRENTSET] = bemobil_interp_avref( EEG , ALLEEG, CURRENTSET, channels_to_interpolate, out_filename, out_filepath)
%
% Inputs:
%   ALLEEG                  - complete EEGLAB data set structure
%   EEG                     - current EEGLAB EEG structure
%   CURRENTSET              - index of current EEGLAB EEG structure within ALLEEG
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
%   EEGLAB, bemobil_interp_reref, bemobil_interp_avref_copy_spatial_filter, bemobil_copy_spatial_filter, pop_interp, pop_reref, pop_interp, 
% 
% Authors: Lukas Gehrke, Marius Klug, 2021

function [ALLEEG, EEG, CURRENTSET] = bemobil_interp( EEG , ALLEEG, CURRENTSET, channels_to_interpolate, out_filename, out_filepath)

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

% Interpolate channels with spherical interpolation

if isempty(channels_to_interpolate)
    disp('No channel indices provided. Attempting to interpolate missing channels from urchanlocs...');
    if ~isempty(EEG.urchanlocs)
        EEG = pop_interp(EEG, EEG.urchanlocs, 'spherical');
        disp('...done.')
        EEG = eeg_checkset(EEG);
    else
        warning('...no urchanlocs present in dataset. Cannot interpolate.');
    end
else
    disp('Interpolating channels that are indicated...');
    EEG = pop_interp(EEG, channels_to_interpolate, 'spherical');
    disp('...done');
end

EEG.etc.interpolated_channels = channels_to_interpolate;

if ~isfield(EEG.etc,'rank') 
    % rank is the number of channels, since we have full rank averef, minus the number of interpolated channels
    EEG.etc.rank = EEG.nbchan - length(channels_to_interpolate);
else
    % rank is reduced by the number of interpolated channels
    EEG.etc.rank = EEG.etc.rank - length(channels_to_interpolate);
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
