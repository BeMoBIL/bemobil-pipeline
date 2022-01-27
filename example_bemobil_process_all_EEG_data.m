
% initialize EEGLAB 
if ~exist('ALLEEG','var')
	eeglab;
end

% initialize FieldTrip 
ft_defaults

% load configuration 
example_bemobil_config;

% if present, load additional metadata saved in example_bemobil_bids_metadata.m
% this is primarily for enhancing documentation and data sharing 
% (no influence on processing if this step is skipped)
example_bemobil_bids_metadata; 

% enter all subjects to process here (you can split it up in more MATLAB instances if you have more CPU power and RAM)
subjects = [64,66,76,78]; 

% set to 1 if all files should be computed, independently of whether they are present on disk or not
force_recompute = 0; 

%% Import 
% (no looping over subjects - enter the whole array of IDs)

% step 1 : convert .xdf to bids
% bemobil_xdf2bids(bemobil_config, subjects) for minimal use 
bemobil_xdf2bids(bemobil_config, subjects, 'general_metadata', general_info, 'motion_metadata', motion_info, 'eeg_metadata',...
    eeg_info, 'participant_metadata', subject_info)

% step 2 : convert bids to .set
bemobil_bids2set(bemobil_config, subjects);


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
    % make sure the data is stored in double precision, large datafiles are supported, no memory mapped objects are
    % used but data is processed locally, and two files are used for storing sets (.set and .fdt)
	try 
        pop_editoptions('option_saveversion6', 0, 'option_single', 0, 'option_memmapdata', 0, 'option_savetwofiles', 1);
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
	[ALLEEG, EEG_preprocessed, CURRENTSET] = bemobil_process_all_EEG_preprocessing(subject, bemobil_config, ALLEEG, EEG, CURRENTSET, force_recompute);

    % start the processing pipeline for AMICA
	bemobil_process_all_AMICA(ALLEEG, EEG_preprocessed, CURRENTSET, subject, bemobil_config, force_recompute);

end

subjects
subject

disp('PROCESSING DONE! YOU CAN CLOSE THE WINDOW NOW!')
