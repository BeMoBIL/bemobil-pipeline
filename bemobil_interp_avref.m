% bemobil_interp_avref() - Interpolates missing channels with spherical
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

function [ EEG ] = bemobil_interp_avref( EEG )

if nargin < 1
	help bemobil_interp_avref;
	return;
end;

% check if preprocessed file already exist and break if it does
out_filename = ['interp_and_averefed_' EEG.filename];
dir_files = dir(EEG.filepath);
if ismember(out_filename, {dir_files.name})
    disp(['Warning: preprocessed file already exists in: ' EEG.filepath '. ' 'Exiting...']);
    return;
end

% Interpolate channels with spherical interpolation
if ~isempty(EEG.urchanlocs)
    EEG = pop_interp(EEG, EEG.urchanlocs, 'spherical');
    disp('Interpolation done.')
    EEG = eeg_checkset(EEG);
else
    disp('No original chanlocs present in dataset. Cannot interpolate.');
    return;
end

% Compute average reference
EEG = pop_reref( EEG, []);

%save data and stop function so manual channel rejection is possible
EEG = eeg_checkset(EEG);
pop_saveset(EEG, strcat(EEG.filepath, out_filename));

disp('Rereferencing done.');