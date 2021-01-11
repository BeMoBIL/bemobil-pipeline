% bemobil_signal_decomposition() - Computes a spatial filter for the EEG data set, which decomposes
% the data into components (e.g. statistically independent components using AMICA - by default).
% AMICA sometimes crashes bevore the first iteration with some combinations of data set length and
% number of threads. Therefore, if this happens and the according error message of WINDOWS is
% closed, AMICA is automatically restarted with one thread less. In case the number of threads are
% reduced to 0, the standard EEGLAB runica() is started.
%
% Usage:
%   >>  [ALLEEG EEG CURRENTSET] = bemobil_signal_decomposition(ALLEEG, EEG, CURRENTSET, amica, numb_models, maxx_threads, data_rank, other_algorithm)
%   >>  [ALLEEG EEG CURRENTSET] = bemobil_signal_decomposition(ALLEEG, EEG, CURRENTSET, amica, numb_models, maxx_threads, data_rank, other_algorithm, out_filename, out_filepath)
%
% Inputs:
%   ALLEEG                  - complete EEGLAB data set structure
%   EEG                     - current EEGLAB EEG structure
%   CURRENTSET              - index of current EEGLAB EEG structure within ALLEEG
%   amica                   - Boolean value (1/0) to use amica or not
%   num_models              - number of models to learn, default is 1
%   max_threads             - maximum of CPU threads to be used for AMICA
%   data_rank               - rank of the data matrix (= number of channels minus number of interpolated
%       channels minus 1, if average referenced)
%   other_algorithm         - currently (20.6.2017) not yet implemented, hopefully SSD and JD will
%       be added here some day
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
%    EEGLAB, runamica15
%
% Authors: Lukas Gehrke, Marius Klug, 2017

function [ALLEEG EEG CURRENTSET] = bemobil_signal_decomposition(ALLEEG, EEG, CURRENTSET,...
	amica, numb_models, maxx_threads, data_rank, other_algorithm, out_filename, out_filepath, AMICA_autoreject, AMICA_n_rej)

% only save a file on disk if both a name and a path are provided
save_file_on_disk = (exist('out_filename', 'var') && exist('out_filepath', 'var') && ...
	~isempty(out_filename) && ~isempty(out_filepath));

% check if file already exist and show warning if it does
if save_file_on_disk
    mkdir(out_filepath); % make sure that folder exists, nothing happens if so
    dir_files = dir(out_filepath);
    if ismember(out_filename, {dir_files.name})
        warning([out_filename ' file already exists in: ' out_filepath '. File will be overwritten...']);
    end
else
    % necessary for AMICA write-out, there won't be other data savings than the AMICA folder!
    out_filepath = EEG.filepath;
    out_filename = EEG.filename;
end

if amica
	
	if ~exist('AMICA_autoreject', 'var')
		AMICA_autoreject = 0;
	end
	
	if ~exist('AMICA_n_rej', 'var')
		AMICA_n_rej = 3;
	end
	
    if isfield(EEG,'datfile') && length(EEG.datfile) > 0
        disp('Found datfile.');
        data = [EEG.filepath '\' EEG.datfile];
        
        
    else
        disp('No datfile field found in EEG structure. Will write temp file in current directory.');
        data = EEG.data(:,:);
        
    end
    
    % delete potentially preexistent folder since it will interfere in case AMICA crashes
    if ~isempty(dir([out_filepath '\' out_filename '_AMICA'])); rmdir([out_filepath '\' out_filename '_AMICA'],'s'); end
    
    disp('Starting AMICA...');
    while maxx_threads > 0
        % try/catch loop because AMICA can crash dependent on the data set and the number of threads
        try
            [w, s, mods] = runamica15(data,...
                'num_models', numb_models,...
                'max_threads', maxx_threads,...
                'outdir', [out_filepath '\' out_filename '_AMICA'],...
                'num_chans', EEG.nbchan,...
                'writestep', 2000,...
                'pcakeep',data_rank,...
				'do_reject',AMICA_autoreject,...
				'numrej',AMICA_n_rej,...
                'write_nd',1,...
                'do_history',0,...
                'histstep',2,...
                'min_dll',1e-9,...
                'min_grad_norm',0.0000005);
            disp('AMICA successfull, storing weights and sphere.');
            EEG.etc.spatial_filter.algorithm = 'AMICA';
            EEG.etc.spatial_filter.AMICAmods = mods;
            
            % if successful, get out of the loop
            break 
            
        catch
            
            % if error, reduce threads by one
            maxx_threads = maxx_threads - 1;
            warning(['AMICA crashed. Reducing maximum threads to ' num2str(maxx_threads)]);
            
        end
    end

    if maxx_threads == 0
    
        error('AMICA crashed with all possible maximum thread options. Try increasing the maximum usable threads of your CPU. If the maximum number of threads has already been tried, you''re pretty much fucked. Ask Jason Palmer, the creator of AMICA.');
        disp('Continuing with default EEGLAB runica() ...');
        amica = 0;
        other_algorithm = 'runica';
    
    end
    
    
end

if ~amica
    % this can't be as the else statement, because the amica can fail and
    % be 0 after having been 1 before
    
    if isempty(other_algorithm)
        other_algorithm = 'runica';
    end
    
    % do other algorithm decomposition
    switch other_algorithm
        
        case 'runica'
            
            [w,s] = runica(EEG.data);
            disp('runica successfull, storing weights and sphere.');
            EEG.etc.spatial_filter.algorithm = 'RUNICA';
        otherwise
            error('Something''s fucky in the other_algorithm loop of bemobil_signal_decomposition...')
            
    end
    
end


% store the actual weights and sphere values in the EEG data set and calculate the rest of the
% spatial filter stuff
EEG.icaweights = w;
EEG.icasphere = s;
EEG = eeg_checkset(EEG);

EEG.etc.spatial_filter.original_data_path = [out_filepath '\' out_filename];

% new data set in EEGLAB
[ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, CURRENTSET, 'gui', 'off');
EEG = eeg_checkset( EEG );

% save on disk
if save_file_on_disk
    EEG = pop_saveset( EEG, 'filename',out_filename,'filepath', out_filepath);
    disp('...done');
end

[ALLEEG EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);