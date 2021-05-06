% close all; clear

example_bemobil_config;

% enter all subjects to process here (you can split it up in more MATLAB instances if you have more CPU power and RAM)
subjects = [1:40]; 
force_recompute = 0;

%% 

if ~exist('ALLEEG','var')
	eeglab;
end

%% processing loop

for subject = subjects
    
    %% prepare filepaths and check if already done
    
	disp(['Subject #' num2str(subject)]);
    
	STUDY = []; CURRENTSTUDY = 0; ALLEEG = [];  CURRENTSET=[]; EEG=[]; EEG_interp_avref = []; EEG_single_subject_final = [];
	
	input_filepath = [bemobil_config.study_folder bemobil_config.raw_EEGLAB_data_folder bemobil_config.filename_prefix num2str(subject)];
	output_filepath = [bemobil_config.study_folder bemobil_config.single_subject_analysis_folder bemobil_config.filename_prefix num2str(subject)];
	
	try
		% load completely processed file
		EEG_single_subject_final = pop_loadset('filename', [ bemobil_config.filename_prefix num2str(subject)...
			'_' bemobil_config.single_subject_cleaned_ICA_filename], 'filepath', output_filepath);
    catch
        disp('...failed. Computing now.')
	end
	
	if ~force_recompute && exist('EEG_single_subject_final','var') && ~isempty(EEG_single_subject_final)
		clear EEG_single_subject_final
		disp('Subject is completely preprocessed already.')
		continue
    end
	
	%% load data that is provided by the BIDS importer
    
    % make sure the data is stored in double precision, large datafiles are supported, and no memory mapped objects are
    % used but data is processed locally
	try 
        pop_editoptions( 'option_saveversion6', 0, 'option_single', 0, 'option_memmapdata', 0);
    catch
        warning('Could NOT edit EEGLAB memory options!!'); 
    end
    
    % load files that were created from xdf to BIDS to EEGLAB
    EEG = pop_loadset('filename',[ bemobil_config.filename_prefix num2str(subject) '_' bemobil_config.merged_filename],'filepath',input_filepath);
    
    %% individual EEG processing to remove non-exp segments
    % it is stongly recommended to remove these because they may contain strong artifacts that confuse channel detection
    % and AMICA
    
    %% processing wrappers for basic processing and AMICA
    
    % do basic preprocessing, line noise removal, and channel interpolation
	[ALLEEG, EEG_interp_avRef, CURRENTSET] = bemobil_process_all_basic(subject, bemobil_config, ALLEEG, EEG, CURRENTSET, force_recompute);

    % start the processing pipeline for AMICA
	bemobil_process_all_AMICA(ALLEEG, EEG_interp_avRef, CURRENTSET, subject, bemobil_config, force_recompute);

end

subjects
subject

disp('PROCESSING DONE! YOU CAN CLOSE THE WINDOW NOW!')
