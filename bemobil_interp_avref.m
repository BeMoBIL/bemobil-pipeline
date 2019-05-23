% bemobil_interp_avref() - Interpolates missing channels with spherical interpolation and
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
%   EEGLAB, bemobil_interp_avref_copy_spatial_filter, bemobil_copy_spatial_filter, pop_interp, pop_reref, pop_interp, 
% 
% Authors: Lukas Gehrke, Marius Klug, 2017

function [ALLEEG, EEG, CURRENTSET] = bemobil_interp_avref( EEG , ALLEEG, CURRENTSET, channels_to_interpolate, out_filename, out_filepath)

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
    EEG.etc.interpolated_channels = channels_to_interpolate;
end

% Compute average reference

EEG_channels_bool = strcmp({EEG.chanlocs.type},'EEG');
EEG_channels = 1:EEG.nbchan;
EEG_channels = EEG_channels(EEG_channels_bool);


if ~any(EEG_channels_bool)
    warning('No channel types in EEG.chanlocs.type are ''EEG''. Continuing with average reference to all channels!')
    EEG = pop_reref( EEG, []);
	 EEG.etc.bemobil_reref = [];
else
    EEG = pop_reref( EEG, EEG_channels,'keepref','on');
	 EEG.etc.bemobil_reref = EEG_channels;
end
disp('Rereferencing done.');

% new data set in EEGLAB
[ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, CURRENTSET, 'gui', 'off');
EEG = eeg_checkset( EEG );

% save on disk
if save_file_on_disk
    EEG = pop_saveset( EEG, 'filename',out_filename,'filepath', out_filepath);
    disp('...done');
end

[ALLEEG EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);
