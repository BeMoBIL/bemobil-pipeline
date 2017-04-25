% bemobil_interp_avref_copyICA() - Interpolates missing channels with spherical
% interpolation and rereferences data to average reference.
%
% Usage:
%   >>  [ EEG ] = bemobil_interp_avref( EEG )
%
% Inputs:
%   EEG     - EEGLAB EEG structure
%    
% Outputs:
%   EEG     - average referenced and channel interpolated EEGLAB EEG structure
%
% See also: 
%   POP_REREF, POP_INTERP, EEGLAB

function [ALLEEG, EEG, CURRENTSET] = bemobil_interp_avref_copyICA( EEG , ALLEEG, CURRENTSET, EEG_index_to_copy_ICA, channels_to_interpolate, out_filename, out_filepath)

if nargin < 1
	help bemobil_interp_avref;
	return;
end;

if ~exist('out_filename', 'var') out_filename = 'interpolated_avRef_ICA_weights.set'; end;
if ~exist('out_filepath', 'var') out_filepath = EEG.filepath; end;

% make sure output folder exists, nothing changes, if yes
mkdir(out_filepath);

% check if preprocessed file already exist and break if it does
dir_files = dir(out_filepath);
if ismember(out_filename, {dir_files.name})
    error(['Warning: ' out_filename ' file already exists in: ' out_filepath '. ' 'Exiting...']);
    %return; use only if warning is provided only on console with disp
end

% Interpolate channels with spherical interpolation

if isempty(channels_to_interpolate)
    disp('No channel indices provided. Attempting to interpolate missing channels from urchanlocs...');
    if ~isempty(EEG.urchanlocs)
        EEG = pop_interp(EEG, EEG.urchanlocs, 'spherical');
        disp('...done.')
        EEG = eeg_checkset(EEG);
    else
        disp('...no urchanlocs present in dataset. Cannot interpolate.');
        return;
    end
else
    disp('Interpolating channels that are indicated...');
    EEG = pop_interp(EEG, channels_to_interpolate, 'spherical');
    disp('...done');
    EEG.etc.interpolated_channels = channels_to_interpolate;
end

% Compute average reference
EEG = pop_reref( EEG, []);
disp('Rereferencing done.');

% Copy ICA weights
if ~isempty(EEG_index_to_copy_ICA)
    disp('Copying ICA weights from provided data set.');
    EEG.icaweights = ALLEEG(EEG_index_to_copy_ICA).icaweights;
    EEG.icawinv = ALLEEG(EEG_index_to_copy_ICA).icawinv;
    EEG.icasphere = ALLEEG(EEG_index_to_copy_ICA).icasphere;
    EEG.icaact = ALLEEG(EEG_index_to_copy_ICA).icaact;
    EEG.icasplinefile = ALLEEG(EEG_index_to_copy_ICA).icasplinefile;
    EEG.icachansind = ALLEEG(EEG_index_to_copy_ICA).icachansind;
    
    EEG = eeg_checkset( EEG );
else
    disp('No data set to copy ICA weights from is provided. Skipping this step.');
end

% save data set
[ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, CURRENTSET, 'gui', 'off');
EEG = eeg_checkset( EEG );
EEG = pop_saveset( EEG, 'filename',out_filename,'filepath', out_filepath);
disp('...done');
[ALLEEG EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);

