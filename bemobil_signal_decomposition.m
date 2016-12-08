% signal_decomposition() - Prepares data for dipole fitting and runs the dipole fitting procedure

% Usage:
%   >>  [ EEG ] = bemobil_dipfit( EEG, headmodel,...
% channels_to_include, components_to_fit, RV_threshold,...
% remove_outside_head, fit_bilateral_dipoles)
%
% Inputs:
%   EEG                   - eeglab EEG struct
%   iteration             - Number, e.g. 1, of current iteration.
%   amica                 - Boolean value (1/0) to use amica or not
%   num_models            - number of models to learn, default is 1
%   max_threads           - see help runamica15
%   other_algorithm       - currently (06.12.2016) not yet implemented
%    
% Outputs:
%   EEG     - EEGLAB EEG structure with ICA weight and sphere matrix
%
% See also: 
%   POP_BEMOBIL_SIGNAL_DECOMPOSITION, POP_RUNAMICA, RUNAMICA15, EEGLAB

function [ EEG ] = bemobil_signal_decomposition( EEG, iteration, amica, numb_models, maxx_threads, other_algorithm)

if nargin < 1
	help bemobil_signal_decomposition;
	return;
end;

% check if decomposed file already exists and break if it does
out_filename = strcat('amica_result_after_', iteration, '_iterations_', EEG.filename);
amica_outdir = strcat(EEG.filepath, 'amica_iteration_', iteration, '/');

dir_files = dir(EEG.filepath);
if ismember(out_filename, {dir_files.name})
    disp(['Warning: preprocessed file already exists in: ' EEG.filepath '. ' 'Exiting...']);
    return;
end

if amica
    [weights, sphere, mods] = runamica15(EEG.data,...
        'num_models', numb_models,...
        'max_threads', maxx_threads,...
        'outdir', amica_outdir);

    %mods = loadmodout15(amica_outdir);
    EEG.icaweights = weights;
    EEG.icasphere = sphere;
    %EEG.ica_mod_prob = mods.mod_prob;
    %EEG.ica_comp_var_explained = mods.svar;
end

if ~isempty(other_algorithm)
    % do other algorithm decomposition
end

%save data and stop function so manual channel rejection is possible
EEG = eeg_checkset(EEG);
pop_saveset(EEG, strcat(EEG.filepath, out_filename));    
end