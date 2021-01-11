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
end
EEG.etc.interpolated_channels = channels_to_interpolate;

% Compute average reference for all EEG channels

EEG_channels_bool = strcmp({EEG.chanlocs.type},'EEG');
REF_channels_bool = strcmp({EEG.chanlocs.type},'REF');
EEG_channels = 1:EEG.nbchan;
EEG_channels = EEG_channels(EEG_channels_bool | REF_channels_bool);


if ~any(EEG_channels_bool)
	
	EEG_channels = [];
	
end


if isempty(EEG.chanlocs(1).ref)
	% no ref was declared during preprocessing, use full rank averef
	% (without needing dependency) Apply average reference after adding
	% initial reference, see fullrankaveref by Makoto Miakoshi (2017) This
	% adds an empty new channel, rereferences, then removes the excess
	% channel again, so the rank is still intact. The reference channel,
	% however, is gone and can't be used fr analyses, so in case Cz was
	% used as reference during recording, it is inaccessible.
	EEG.nbchan = EEG.nbchan+1;
	EEG.data(end+1,:) = zeros(1, EEG.pnts);
	EEG.chanlocs(1,EEG.nbchan).labels = 'initialReference';
	EEG = pop_reref( EEG, EEG_channels,'keepref','on');
	EEG = pop_select( EEG,'nochannel',{'initialReference'});

    % rank is the number of channels, since we have full rank averef, minus the number of interpolated channels
    EEG.etc.rank = EEG.nbchan - length(channels_to_interpolate);
else
	% ref was declared, keep it as channel. this means we have an extra
	% channel, e.g. 129 instead of 128 electrodes, and the former reference
	% carries information. however, the original rank is still EEG.nbchan - 1, so 128!
	EEG = pop_reref( EEG, EEG_channels,'keepref','on');
    
    % rank is the number of channels - 1, minus the number of interpolated channels
    EEG.etc.rank = EEG.nbchan - 1 - length(channels_to_interpolate);
end

EEG.etc.bemobil_reref = EEG_channels;
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
