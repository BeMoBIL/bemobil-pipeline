% bemobil_signal_decomposition() - Computes a spatial filter for the EEG data set, which decomposes the data into
% components (e.g. statistically independent components using AMICA - by default). AMICA sometimes crashes before the
% first iteration with some combinations of data set length and number of threads. Therefore, if this happens AMICA is
% automatically restarted with one thread less. Allows automatic rejection of samples in AMICA, which effectively
% replaces time-domain artifact cleaning.
%
% Usage:
%   >>  [ALLEEG EEG CURRENTSET] = bemobil_signal_decomposition(ALLEEG, EEG, CURRENTSET,...
%           amica, numb_models, maxx_threads, data_rank, other_algorithm, out_filename, out_filepath, AMICA_autoreject,...
%           AMICA_n_rej, AMICA_reject_sigma_threshold, AMICA_max_iter)

% Inputs:
%   ALLEEG                          - complete EEGLAB data set structure
%   EEG                             - current EEGLAB EEG structure
%   CURRENTSET                      - index of current EEGLAB EEG structure within ALLEEG
%   amica                           - Boolean value (1/0) to use amica or not
%   numb_models                     - number of models to learn, default is 1
%   maxx_threads                    - maximum of CPU threads to be used for AMICA
%   data_rank                       - rank of the data matrix (= number of channels minus number of interpolated
%       channels minus 1, if average referenced)
%   other_algorithm                 - currently (20.6.2017) not yet implemented, hopefully SSD and JD will
%       be added here some day
%   out_filename                    - output filename (OPTIONAL)
%   out_filepath                    - output filepath (OPTIONAL - File will only be saved on disk
%       if both a name and a path are provided)
%   AMICA_autoreject                - flag for doing rejection of time points, def=1 (OPTIONAL)
%   AMICA_n_rej                     - for rejection, number of rejections to perform, def=5 (OPTIONAL)
%   AMICA_reject_sigma_threshold    - for rejection, sigma threshold of log likelyhood of samples to reject, def=3 (OPTIONAL)
%   AMICA_max_iter                  - maximum number of iterations for AMICA
%
% Outputs:
%   ALLEEG                          - complete EEGLAB data set structure
%   EEG                             - current EEGLAB EEG structure
%   Currentset                      - index of current EEGLAB EEG structure within ALLEEG
%
%   .set data file of current EEGLAB EEG structure stored on disk (OPTIONALLY)
%
% See also:
%    EEGLAB, runamica15
%
% Authors: Lukas Gehrke, Marius Klug, 2021

function [ALLEEG EEG CURRENTSET] = bemobil_signal_decomposition(ALLEEG, EEG, CURRENTSET,...
	amica, numb_models, maxx_threads, data_rank, other_algorithm, out_filename, out_filepath, AMICA_autoreject,...
    AMICA_n_rej, AMICA_reject_sigma_threshold, AMICA_max_iter)

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

amica_crashed = false; 

if amica
	
	if ~exist('AMICA_autoreject', 'var')
		AMICA_autoreject = 1;
	end
	
    if ~exist('AMICA_n_rej', 'var')
		AMICA_n_rej = 5;
    end
    
    if ~exist('AMICA_reject_sigma_threshold', 'var')
		AMICA_reject_sigma_threshold = 3;
    end
    
    if ~exist('AMICA_max_iter', 'var')
		AMICA_max_iter = 2000;
    end
	
    if isfield(EEG,'datfile') && ~isempty(EEG.datfile)
        disp('Found datfile.');
        data = fullfile(EEG.filepath, EEG.datfile);
        
    else
        disp('No datfile field found in EEG structure. Will write temp file in current directory.');
        data = EEG.data(:,:);
        
    end
    
    % delete potentially preexistent folder since it will interfere in case AMICA crashes
    try
        rmdir(fullfile(out_filepath, [erase(out_filename,'.set') '_outfiles']),'s'); 
    end
    
    disp('Starting AMICA...');
    while maxx_threads > 0
        
        % try/catch loop because AMICA can crash dependent on the data set and the number of threads
        try
            [w, s, mods] = runamica15(data,...
                'num_models', numb_models,...
                'max_threads', maxx_threads,...
                'outdir', fullfile(out_filepath, [erase(out_filename,'.set') '_outfiles']),...
                'num_chans', EEG.nbchan,...
                'max_iter',AMICA_max_iter,...
                'writestep', AMICA_max_iter,...
                'pcakeep',data_rank,...
				'do_reject',AMICA_autoreject,...
				'numrej',AMICA_n_rej,...
                'rejsig',AMICA_reject_sigma_threshold,...
                'write_nd',1,...
                'do_history',0,...
                'histstep',2,...
                'min_dll',0.000000001,...
                'min_grad_norm',0.0000005);
            disp('AMICA successfull, storing weights and sphere.');
            EEG.etc.spatial_filter.algorithm = 'AMICA';
            EEG.etc.spatial_filter.AMICAmods = mods;
            
            % add information about AMICA autorejected time points
            sample_mask = EEG.etc.spatial_filter.AMICAmods.Lt == 0;
            EEG.etc.bad_samples = sample_mask;
            EEG.etc.bad_samples_percent = sum(EEG.etc.bad_samples) / length(EEG.etc.bad_samples) * 100;

            % find latency of regions
            remove_data_intervals = reshape(find(diff([false sample_mask false])),2,[])';
            remove_data_intervals(:,2) = remove_data_intervals(:,2)-1;
            EEG.etc.remove_data_intervals = remove_data_intervals;
            
            % if successful, get out of the loop
            break 
            
        catch error_message
            
            % if error, reduce threads by one
            maxx_threads = maxx_threads - 1;
            warning(['AMICA crashed. This could be wrong input or a bug in AMICA depending on the number of threads. Reducing maximum threads to ' num2str(maxx_threads)]);
            
        end
    end

    if maxx_threads == 0
    
        % throw error message from runamica15 as warning
        warning(error_message.message)
        
        % potential issue with number of threads
        error('AMICA crashed with all possible threading options. Is AMICA correctly installed? Does the filepath contain NO spaces? Are all inputs to the function correct? If yes, please contact the BeMoBIL.');        
        
%         disp('Continuing with default EEGLAB runica() ...');
%         amica_crashed   = true; 
%         other_algorithm = 'runica';
    
    end
    
end

if ~amica || amica_crashed

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
            error('Invalid definition of input other_algorithm.')
            
    end
    
end


% store the actual weights and sphere values in the EEG data set and calculate the rest of the
% spatial filter stuff
EEG.icaweights = w;
EEG.icasphere = s;
EEG = eeg_checkset(EEG);

EEG.etc.spatial_filter.original_data_path = fullfile(out_filepath, out_filename);

% new data set in EEGLAB
[ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, CURRENTSET, 'gui', 'off');
EEG = eeg_checkset( EEG );

% save on disk
if save_file_on_disk
    EEG = pop_saveset( EEG, 'filename',out_filename,'filepath', out_filepath);
    disp('...done');
end

[ALLEEG EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);