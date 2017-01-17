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

function [ALLEEG EEG CURRENTSET] = bemobil_signal_decomposition(ALLEEG, EEG, CURRENTSET, iteration, amica, numb_models, maxx_threads, other_algorithm)

if nargin < 1
    help bemobil_signal_decomposition;
    return;
end;

% check if decomposed file already exists and break if it does
out_filename = ['postICA' num2str(iteration)];
amica_outdir = strcat(EEG.filepath, 'amica_iteration_', iteration);

dir_files = dir(EEG.filepath);
if ismember(out_filename, {dir_files.name})
    disp(['Warning: preprocessed file already exists in: ' EEG.filepath '. ' 'Exiting...']);
    return;
end

if amica
    if isfield(EEG,'datfile') && length(EEG.datfile) > 0
        disp('Found datfile');
        while maxx_threads > 0
            try
                [w, s, mods] = runamica15([EEG.filepath EEG.datfile],...
                    'num_models', numb_models,...
                    'max_threads', maxx_threads,...
                    'outdir', amica_outdir,...
                    'num_chans', EEG.nbchan);
                break
            catch
                maxx_threads = maxx_threads - 1;
                warning(['AMICA crashed. Reducing maximum threads to ' num2str(maxx_threads)]);
                
            end
        end
        warning('AMICA crashed with all possible maximum thread options. Try increasing the maximum usable threads. If the maximum number of threads has already been tried, you''re pretty much fucked. Ask Jason Palmer, the creator of AMICA.');
        disp('Continuing with default runica...');
        other_algorithm = 'runica';
    else
        disp('No datfile field found in EEG structure. Will write temp file in current directory.');
        while maxx_threads > 0
            try
                [w, s, mods] = runamica15(EEG.data(:,:),...
                    'num_models', numb_models,...
                    'max_threads', maxx_threads,...
                    'outdir', amica_outdir);
                break
            catch
                maxx_threads = maxx_threads - 1;
                warning(['AMICA crashed. Reducing maximum threads to ' num2str(maxx_threads)]);
                
            end
        end
        warning('AMICA crashed with all possible maximum thread options. Try increasing the maximum usable threads. If the maximum number of threads has already been tried, you''re pretty much fucked. Ask Jason Palmer, the creator of AMICA.');
        disp('Continuing with default runica...');
        other_algorithm = 'runica';
    end
    
    
    %mods = loadmodout15(amica_outdir);
    %     EEG.icaweights = w;
    %     EEG.icasphere = s;
    %EEG.ica_mod_prob = mods.mod_prob;
    %EEG.ica_comp_var_explained = mods.svar;
else
    
    other_algorithm = 'runica';
    
end

if ~isempty(other_algorithm)
    % do other algorithm decomposition
    
    if strcmp(other_algorithm, 'runica')
        
        [w,s] = runica(EEG.data);
        %         EEG.icaweights = w;
        %         EEG.icasphere = s;
        
    end
end

EEG.icaweights = w;
EEG.icasphere = s;
EEG = eeg_checkset(EEG);

% saving data set
[ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, CURRENTSET, 'gui', 'off');
EEG = eeg_checkset( EEG );
EEG = pop_saveset( EEG, 'filename',out_filename,'filepath', [ ALLEEG(CURRENTSET-1).filepath '\']);
disp('...done');
[ALLEEG EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);

end