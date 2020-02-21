% bemobil_process_all_AMICA - wrapper function that incorporates all necessary processing steps from the basic EEG
% struct (e.g. all blocks merged together, nothing else done before) up to the finished dataset which has channels
% interpolated and all AMICA information copied. The AMICA is computed on a dataset that made use of automatic channel
% and time domain cleaning. Additionally, information of dipole fitting and automatic IC classification with ICLabel is
% present. A processing config struct is necessary. For an example please see the EEG_processing_example script!
%
% The intermediate files are stored on the disk.
%
% Usage:
%   >>  [ALLEEG, EEG_single_subject_final, CURRENTSET] = bemobil_process_all_AMICA(ALLEEG, EEG_to_process, CURRENTSET, subject, bemobil_config)
%
% Inputs:
%   ALLEEG                    - complete EEGLAB data set structure
%   EEG_to_process            - EEGLAB EEG structure that should be processed. Best to have all blocks merged into one
%                                file.
%   CURRENTSET                - index of current EEGLAB EEG structure within ALLEEG
%   subject                   - subject number of the current subject (necessary for filepaths and storage)
%   bemobil_config            - configuration struct with all necessary information. See EEG_processing_example file
%                                that comes with this function!
%
% Outputs:
%   ALLEEG                    - complete EEGLAB data set structure
%   EEG_single_subject_final  - current EEGLAB EEG structure
%   Currentset                - index of current EEGLAB EEG structure within ALLEEG
%
%   .set data file of current EEGLAB EEG structure stored on disk (OPTIONALLY)
%
% See also:
%   EEGLAB
%
% Authors: Marius Klug, 2019

function [ALLEEG, EEG_single_subject_final, CURRENTSET] = bemobil_process_all_AMICA(ALLEEG, EEG_interp_avRef, CURRENTSET, subject, bemobil_config, force_recompute)

if ~exist('force_recompute','var')
	force_recompute = false;
end
if force_recompute
	warning('RECOMPUTING OLD FILES IF FOUND!!!')
end

% check if the entire processing was done already
output_filepath = [bemobil_config.study_folder bemobil_config.single_subject_analysis_folder bemobil_config.filename_prefix num2str(subject)];
	
if ~force_recompute
	try
		
		EEG = pop_loadset('filename', [bemobil_config.filename_prefix num2str(subject) '_'...
			bemobil_config.copy_weights_interpolate_avRef_filename], 'filepath', output_filepath);
		[ALLEEG, EEG_single_subject_final, CURRENTSET] = pop_newset(ALLEEG, EEG, 0,'study',0);
		
		warning('Old cleaned AMICA file already existed, using that file!')
		
	end
end

if ~exist('EEG_single_subject_final','var')
	
	%% highpass filter for AMICA
	
	output_filepath = [bemobil_config.study_folder bemobil_config.spatial_filters_folder...
		bemobil_config.spatial_filters_folder_AMICA bemobil_config.filename_prefix num2str(subject)];
	
	% check if the part was done already
	if ~force_recompute
		try
			
			EEG = pop_loadset('filename', [bemobil_config.filename_prefix num2str(subject) '_'...
				bemobil_config.filtered_filename], 'filepath', output_filepath);
			[ALLEEG, EEG_filtered_for_AMICA, CURRENTSET] = pop_newset(ALLEEG, EEG, 0,'study',0);
			
			warning('Old filtered file already existed, using that file!')
			
		end
	end
	
	if ~exist('EEG_filtered_for_AMICA','var')
		
		disp('Filtering the data for AMICA...')
		
		% delete events to save disk space and RAM 
% 		EEG_interp_avRef.event = []; DONT DO THIS, IT'S THE DATA SET THAT GET'S COPIED LATER!
		
		% filters the data set separately for low and high cutoff frequencies, stores all relevant info in the set
		[ALLEEG, EEG_filtered_for_AMICA, CURRENTSET] = bemobil_filter(ALLEEG, EEG_interp_avRef, CURRENTSET,...
			bemobil_config.filter_lowCutoffFreqAMICA, bemobil_config.filter_highCutoffFreqAMICA,...
			[bemobil_config.filename_prefix num2str(subject) '_' bemobil_config.filtered_filename], output_filepath);
		
	end
	
	%% AMICA Loop 1
	% this is only to find eye components that should be
	% projected out before automatic continuous data cleaning
	
	output_path = [bemobil_config.study_folder bemobil_config.spatial_filters_folder...
		bemobil_config.spatial_filters_folder_AMICA];
	
	output_filepath = [output_path bemobil_config.filename_prefix num2str(subject)];
	
	% check if the part was done already
	if ~force_recompute
		try
			
			EEG = pop_loadset('filename', [bemobil_config.filename_prefix num2str(subject) '_'...
				bemobil_config.amica_raw_filename_output], 'filepath', output_filepath);
			[ALLEEG, EEG_AMICA_raw, CURRENTSET] = pop_newset(ALLEEG, EEG, 0,'study',0);
			
			warning('Old raw AMICA file already existed, using that file!')
			
		end
	end
	
	if ~exist('EEG_AMICA_raw','var')
		
		% rank of the data set is reduced by 1, because of average reference,
		% and then reduced by the number of interpolated channels
		if ~isfield(EEG_filtered_for_AMICA.etc, 'interpolated_channels')
			warning('No channels were interpolated. Is this correct?');
			rank = EEG_filtered_for_AMICA.nbchan - 1;
		else
			rank = EEG_filtered_for_AMICA.nbchan - 1 - length(EEG_filtered_for_AMICA.etc.interpolated_channels);
		end
		
		% running signal decomposition with values specified above
		disp('AMICA computation on raw data for projecting eye components out...');
		[ALLEEG, EEG_AMICA_raw, CURRENTSET] = bemobil_signal_decomposition(ALLEEG, EEG_filtered_for_AMICA, ...
			CURRENTSET, true, bemobil_config.num_models, bemobil_config.max_threads, rank, [], ...
			[bemobil_config.filename_prefix num2str(subject) '_' bemobil_config.amica_raw_filename_output], output_filepath);
		
	end
	
	%% Determine (ICLabel) eye components and project out
	
	output_path = [bemobil_config.study_folder bemobil_config.spatial_filters_folder...
		bemobil_config.spatial_filters_folder_AMICA];
	
	output_filepath = [output_path bemobil_config.filename_prefix num2str(subject)];
	
	% check if the part was done already
	if ~force_recompute
		try
			
			EEG = pop_loadset('filename', [bemobil_config.filename_prefix num2str(subject) '_'...
				bemobil_config.amica_chan_no_eye_filename_output], 'filepath', output_filepath);
			[ALLEEG, EEG_AMICA_no_eyes, CURRENTSET] = pop_newset(ALLEEG, EEG, 0,'study',0);
			
			warning('Old no eyes AMICA file already existed, using that file!')
			
		end
	end
	
	if ~exist('EEG_AMICA_no_eyes','var')
		
		% iclabel
		disp('Determine (ICLabel) eye components and project out...');
		EEG_AMICA_raw = iclabel(EEG_AMICA_raw,'lite');
		
		% define eyes
		eyes = find(EEG_AMICA_raw.etc.ic_classification.ICLabel.classifications(:,3) > bemobil_config.eye_threshold);
		
		% plot eyes and save
		h = bemobil_plot_patterns(EEG_AMICA_raw.icawinv,EEG_AMICA_raw.chanlocs,'weights',...
			EEG_AMICA_raw.etc.ic_classification.ICLabel.classifications(:,3),'minweight',bemobil_config.eye_threshold);
		
		savefig(h, [output_filepath '\' bemobil_config.filename_prefix num2str(subject) '_eye-comps']);
		close(gcf);
		
		% project out eyes and save
		EEG_AMICA_no_eyes = pop_subcomp(EEG_AMICA_raw, eyes);
		
		% save RAM
		clear EEG_AMICA_raw
		
		pop_saveset(EEG_AMICA_no_eyes, 'filename', [bemobil_config.filename_prefix num2str(subject) '_'...
			bemobil_config.amica_chan_no_eye_filename_output], 'filepath', output_filepath);
		disp('... done!')
		
	end
	
	if ~isfield(EEG_AMICA_no_eyes.etc,'auto_continuous_cleaning')
		
		% determine automatic time domain cleaning boundaries on the channel level
		
		input_filepath = [bemobil_config.study_folder bemobil_config.spatial_filters_folder...
			bemobil_config.spatial_filters_folder_AMICA bemobil_config.filename_prefix num2str(subject)];
		output_filepath = [bemobil_config.study_folder bemobil_config.spatial_filters_folder...
			bemobil_config.spatial_filters_folder_AMICA bemobil_config.filename_prefix num2str(subject)];
		
		EEG_AMICA_no_eyes.event = [];
		EEG_AMICA_no_eyes = eeg_checkset( EEG_AMICA_no_eyes );
		
		%%% specify folder names
		datapath_specifications.datapath_original_files=input_filepath; %%% keep last \; single subject folder
		datapath_specifications.datapath_save_files=output_filepath;     %%% keep last \; path for saving updated EEG
		datapath_specifications.datapath_save_figures=output_filepath;   %%% keep last \; path for saving figures of cleaning
		
		%%% specify file names
		filename_specifications.file_name_original_EEG=...
			[bemobil_config.filename_prefix num2str(subject) '_' bemobil_config.amica_chan_no_eye_filename_output];   %%% loads "fresh" EEG (raw, unfiltered)
		filename_specifications.filename_saveBadEpochIndices='';
		
		automatic_cleaning_settings.cleaned_data_type='sensor data'; %%% ICA not implemented yet; usually bad segments found on sensor level are also fine for IC later on
		
		%%% select channels that should be considered for cleaning
		automatic_cleaning_settings.selected_sensor_channels_for_cleaning=[]; %%% [] use all available channels for cleaning, else specify [1 2 ...]; currently same channels for all subjects alike
		%%% select channel(s) for cleaning before vs. after
		automatic_cleaning_settings.chan_select_plot_before_after=[5];  %%% [] use all available channels for cleaning
		
		%%% define frequency band (band-pass filter) only for cleaning
		automatic_cleaning_settings.band_artifact_cleaning=[1 40];  %%% in Hz; [] for no filter
		automatic_cleaning_settings.band_stop_artifact_cleaning=[]; %%% [] for nothing; else e.g. [48 52] for removal of line artifacts (notch filter)
		automatic_cleaning_settings.band_filtorder=2;               %%% for IIR Butterworth filter; since filtfilt is used, the effective order is double (here 4)
		
		%%% further settings for the cleaning algorithm
		automatic_cleaning_settings.analyzed_files_N=1;  %%% so far: analyze only 1 file! (no appended different conditions)
		automatic_cleaning_settings.crit_all=bemobil_config.automatic_cleaning_threshold_to_keep; %%% e.g., 0.9=90% keep amount of epochs; 1 value if only 1 file (no appended recordings); else indicate, e.g., 4 values for 4 appended files
		automatic_cleaning_settings.wind_ms=1000;    %%% in ms, epochs for finding large artifacts
		automatic_cleaning_settings.crit_keep_sec=[]; %%% in seconds; for NOT removing any additional data, set: automatic_cleaning_settings.crit_keep_sec=automatic_cleaning_settings.wind_ms/1000; else: value should be multiple of "wind_ms" for additionally removing data, i.e., keep uninterrupted "good" data segments not shorter than this value
		automatic_cleaning_settings.crit_percent_sample_epoch=0.2;  %%% [] for nothing; e.g., 0.1 = 10%; remove epochs if they have more than x % of samples with zero or NaN ("flat line")
		automatic_cleaning_settings.weighting_factor_epoch_cleaning_methods=[1 1 2];     %%% method I mean of epochs, method II = channel heterogeneity --> SD across channels of mean_epochs; method III = channel heterogeneity --> Mahal. distance of mean_epochs across channels; recommended: put method I not at zero, because Mahal. might not yield results if data set is extremely short
		automatic_cleaning_settings.visual_inspection_mode=false;  %%% =false if visual threshold rejection after automatic cleaning should not be applied; in this case, bad segments from previous automatic artifact rejection are taken
		if ~automatic_cleaning_settings.visual_inspection_mode
			automatic_cleaning_settings.threshold_visual_reject=zeros(1,automatic_cleaning_settings.analyzed_files_N);
		end
		
		% 2-5) loading, cleaning, and saving (results, figures)
		%%% only indices are saved in the struct "auto_continuous_cleaning"
		%%% so far: before/after comparisons are saved as .jpg only, because
		%%% otherwise usually files are too large with continuous data. However, if
		%%% you open the "wrapper_automatic_cleaning_continuous_EEG" and run each
		%%% section manually, you can save the same figure as .fig (it's just disabled by
		%%% now). Please do not enable it in the wrapper itself!
		
		addpath(genpath('P:\Project_Friederike\2017_spot_rotation\1_analysis\analysis_Matlab_diaries_scripts_releases\3_release_internal_use_only'))
		
		disp('Determining continuous data cleaning boundaries...')
		[auto_continuous_cleaning]=wrapper_automatic_cleaning_continuous_EEG(datapath_specifications,filename_specifications,...
			automatic_cleaning_settings);
		
		% copy cleaning results and save dataset
		EEG_AMICA_no_eyes.etc.auto_continuous_cleaning = auto_continuous_cleaning;
		
		% add buffers to bad epochs
		addpath('P:\Lukas_Gehrke\toolboxes\utils_LG');
		
		disp('Adding buffers to data cleaning boundaries...')
		EEG_AMICA_no_eyes = add_leading_trailing_samples_FH_cleaning(EEG_AMICA_no_eyes, bemobil_config.buffer_length);
		
		pop_saveset(EEG_AMICA_no_eyes, 'filename', [bemobil_config.filename_prefix num2str(subject) '_'...
			bemobil_config.amica_chan_no_eye_filename_output], 'filepath', output_filepath);
		disp('... done!')
		close(gcf);
		
	end
	
	% save RAM
	clear EEG_AMICA_raw
	
	%% ICA loop 2
	
	output_filepath = [bemobil_config.study_folder bemobil_config.spatial_filters_folder...
		bemobil_config.spatial_filters_folder_AMICA bemobil_config.filename_prefix num2str(subject)];
	
	% check if the part was done already
	if ~force_recompute
		try
			
			EEG = pop_loadset('filename', [bemobil_config.filename_prefix num2str(subject) '_'...
				bemobil_config.amica_filename_output], 'filepath', output_filepath);
			[ALLEEG, EEG_AMICA_cleaned, CURRENTSET] = pop_newset(ALLEEG, EEG, 0,'study',0);
			
			warning('Old cleaned AMICA file already existed, using that file!')
			
		end
	end
	
	if ~exist('EEG_AMICA_cleaned','var')
		
		% reject data from previously computed FH channel data rejection
		% get cleaning indices
		invalid_segments_index = EEG_AMICA_no_eyes.etc.auto_continuous_cleaning.cleaning_buffer.invalid_segments_final_start_stop_sample_buffered;
		
		% save RAM
		clear EEG_AMICA_no_eyes
		
		% apply rejection on original dataset EEG_interp_avRef
		EEG_cleaned_for_AMICA = eeg_eegrej(EEG_filtered_for_AMICA, invalid_segments_index);
		
		% rank of the data set is reduced by 1, because of average reference,
		% and then reduced by the number of interpolated channels
		if ~isfield(EEG_filtered_for_AMICA.etc, 'interpolated_channels')
			warndlg(['No channels were interpolated, subject = ' num2str(subject) '. Is this correct?']);
			rank = EEG_filtered_for_AMICA.nbchan - 1;
		else
			rank = EEG_filtered_for_AMICA.nbchan - 1 - length(EEG_filtered_for_AMICA.etc.interpolated_channels);
		end
		
		clear EEG_filtered_for_AMICA
		
		% running signal decomposition with values specified above
		
		disp('Final AMICA computation on cleaned data...');
		[ALLEEG, EEG_AMICA_cleaned, CURRENTSET] = bemobil_signal_decomposition(ALLEEG, EEG_cleaned_for_AMICA, ...
			CURRENTSET, true, bemobil_config.num_models, bemobil_config.max_threads, rank, [], ...
			[bemobil_config.filename_prefix num2str(subject) '_' bemobil_config.amica_filename_output], output_filepath);
		
	end
	
	% save RAM
	clear EEG_AMICA_no_eyes
	clear EEG_filtered_for_AMICA
	clear EEG_cleaned_for_AMICA
	
	%% Warping of locations and dipole fitting, plus runing ICLabel
	% renames the specified channels, warps the chanlocs on a standard head model and fits dipoles for
	% each IC below the threshold of residual variance
	
	output_filepath = [bemobil_config.study_folder bemobil_config.spatial_filters_folder...
		bemobil_config.spatial_filters_folder_AMICA bemobil_config.filename_prefix num2str(subject)];
	
	% check if the part was done already
	if ~force_recompute
		try
			
			EEG = pop_loadset('filename', [bemobil_config.filename_prefix num2str(subject) '_'...
				bemobil_config.warped_dipfitted_filename], 'filepath', output_filepath);
			[ALLEEG, EEG_AMICA_final, CURRENTSET] = pop_newset(ALLEEG, EEG, 0,'study',0);
			
			warning('Old cleaned AMICA file already existed, using that file!')
			
		end
	end
	
	if ~exist('EEG_AMICA_final','var')
		
		% compute iclabel scores
		disp('ICLabel component classification...');
		EEG_AMICA_cleaned = iclabel(EEG_AMICA_cleaned,'lite');
		
		% do the warp and dipfit
		disp('Dipole fitting...');
		[ALLEEG, EEG_AMICA_final, CURRENTSET] = bemobil_dipfit( EEG_AMICA_cleaned , ALLEEG, CURRENTSET, bemobil_config.warping_channel_names,...
			bemobil_config.residualVariance_threshold,...
			bemobil_config.do_remove_outside_head, bemobil_config.number_of_dipoles,...
			[bemobil_config.filename_prefix num2str(subject) '_' bemobil_config.warped_dipfitted_filename], output_filepath);
		
	end
	
	% save RAM
	clear EEG_AMICA_cleaned
	
	%% Final step: copy the spatial filter data into the raw full data set for further single subject processing
	
	output_filepath = [bemobil_config.study_folder bemobil_config.single_subject_analysis_folder bemobil_config.filename_prefix num2str(subject)];
	
	disp('Copying all information into full length dataset for single subject processing...');
	[ALLEEG, EEG_single_subject_final, CURRENTSET] = bemobil_copy_spatial_filter(EEG_interp_avRef, ALLEEG, CURRENTSET,...
		EEG_AMICA_final, [bemobil_config.filename_prefix num2str(subject) '_'...
		bemobil_config.copy_weights_interpolate_avRef_filename], output_filepath);
	
end

disp('Entire processing done!');
