% close all; clear

% load configuration 
example_bemobil_config;

% if present, load additional metadata saved in example_bemobil_bids_metadata.m
% this is primarily for enhancing documentation and data sharing 
% (no influence on processing if this step is skipped)
example_bemobil_bids_metadata; 

% enter all subjects to process here (you can split it up in more MATLAB instances if you have more CPU power and RAM)
subjects = [64,66,76,78]; 
force_recomp = false;


%% Import 
% (no looping over subjects - enter the whole array of IDs)

% step 1 : convert .xdf to bids
% bemobil_xdf2bids(bemobil_config, subjects) for minimal use 
bemobil_xdf2bids(bemobil_config, subjects, 'general_metadata', general_info, 'motion_metadata', motion_info, 'eeg_metadata', eeg_info, 'participant_metadata', subject_info)

% step 2 : convert bids to .set
bemobil_set2bids(bemobil_config, subject);


%% processing loop

if ~exist('ALLEEG','var')
	eeglab;
end

if ~exist('mobilab','var') 
	runmobilab;
end
%%

for subject = subjects
    
	STUDY = []; CURRENTSTUDY = 0; ALLEEG = [];  CURRENTSET=[]; EEG=[]; EEG_interp_avref = [];
	
	input_filepath = [bemobil_config.study_folder bemobil_config.single_subject_analysis_folder bemobil_config.filename_prefix num2str(subject)];
	output_filepath = [bemobil_config.study_folder bemobil_config.single_subject_analysis_folder bemobil_config.filename_prefix num2str(subject)];
	
	try
		% load completely processed file
		EEG_single_subject_final = pop_loadset('filename', [ bemobil_config.filename_prefix num2str(subject)...
			'_' bemobil_config.single_subject_cleaned_ICA_filename], 'filepath', output_filepath);
	end
	
	if ~force_recomp && exist('EEG_single_subject_final','var')
		clear EEG_single_subject_final
		disp('Subject is completely preprocessed already.')
		continue
    end
	
	disp(['Subject #' num2str(subject)]);
	
% 	load xdf files and process them with mobilab, export to eeglab, split MoBI and merge all conditions for EEG, do
% 	basic preprocessing and interpolation
	[ALLEEG, EEG_interp_avref, CURRENTSET] = bemobil_process_all_mobilab(subject, bemobil_config, ALLEEG, CURRENTSET, mobilab, force_recomp);

    % get rid of memory mapped object storage by mobilab
	try 
        pop_editoptions( 'option_saveversion6', 0, 'option_single', 0, 'option_memmapdata', 0);
    catch
        warning('Could NOT edit EEGLAB memory options!!'); 
    end
    
% 	start the processing pipeline for AMICA
	bemobil_process_all_AMICA(ALLEEG, EEG_interp_avref, CURRENTSET, subject, bemobil_config, force_recomp);

end

subjects
subject

close all
clear
disp('PROCESSING DONE! YOU CAN CLOSE THE WINDOW NOW!')
