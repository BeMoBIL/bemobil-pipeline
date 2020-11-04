close all; clear

%% ONLY CHANCE THESE PARTS!

subjects = 1;

process_config.study_folder = 'P:\situation_awareness\data\third_real_pilot\';

process_config.filenames = {'block_calibration','block_slow_narrow','block_slow_wide','block_fast_narrow',...
   'block_fast_wide','block_fast_narrow_grey','block_fast_wide_grey','block_fast_narrow_stable_head',...
   'block_fast_narrow_stable_both'};

% enter channels that you did not use at all (e.g. with the MoBI 160 chan layout, only 157 chans are used):
% process_config.channels_to_remove = {'N29' 'N30' 'N31'};
process_config.channels_to_remove = [];

% enter EOG channel names here:
% process_config.eog_channels  = {'G16' 'G32'};
process_config.eog_channels  = {'EOG2'};

% leave this empty if you have standard channel names that should use standard locations:
% process_config.channel_locations_filename = 'channel_locations.elc';
process_config.channel_locations_filename = [];


%% everything from here is according to the general pipeline, changes only recommended if you know the whole structure

% general foldernames and filenames
process_config.raw_data_folder = '0_raw-data\';
process_config.raw_EEGLAB_data_folder = '2_basic-EEGLAB\';
process_config.spatial_filters_folder = '3_spatial-filters\';
process_config.spatial_filters_folder_AMICA = '3-1_AMICA\';
process_config.single_subject_analysis_folder = '4_single_subject_analysis\';

process_config.merged_filename = 'merged.set';
process_config.preprocessed_filename = 'preprocessed.set';
process_config.interpolated_avRef_filename = 'interpolated_avRef.set';
process_config.filtered_filename = 'filtered.set';
process_config.amica_raw_filename_output = 'postAMICA_raw.set';
process_config.amica_chan_no_eye_filename_output = 'preAMICA_no_eyes.set';
process_config.amica_filename_output = 'postAMICA_cleaned.set';
process_config.warped_dipfitted_filename = 'warped_dipfitted.set';
process_config.copy_weights_interpolate_avRef_filename = 'interp_avRef_ICA.set';

% preprocessing
process_config.resample_freq = 250;

%%% AMICA

% on some PCs AMICA may crash before the first iteration if the number of
% threads and the amount the data does not suit the algorithm. Jason Palmer
% has been informed, but no fix so far. just roll with it. if you see the
% first iteration working there won't be any further crashes. in this case
% just press "close program" or the like and the bemobil_spatial_filter
% algorithm will AUTOMATICALLY reduce the number of threads and start AMICA
% again. this way you will always have the maximum number
% of threads that should be used for AMICA. check in the
% task manager how many threads you have theoretically available and think
% how much computing power you want to devote for AMICA. on the bpn-s1
% server, 12 is half of the capacity and can be used. be sure to check with
% either Ole or your supervisor and also check the CPU usage in the task
% manager before!

% 4 threads are most effective for single subject speed, more threads don't
% really shorten the calculation time much. best efficiency is using just 1
% thread and have as many matlab instances open as possible (limited by the
% CPU usage). Remember your RAM limit in this case.


process_config.filter_lowCutoffFreqAMICA = 1;
process_config.filter_highCutoffFreqAMICA = [];
process_config.max_threads = 4;
process_config.num_models = 1;

% warp electrodemontage and run dipfit
process_config.warping_channel_names = [];
process_config.residualVariance_threshold = 100;
process_config.do_remove_outside_head = 'off';
process_config.number_of_dipoles = 1;

% IC_label
process_config.eye_threshold = 0.7;

% FHs cleaning
process_config.buffer_length = 0.49;
process_config.automatic_cleaning_threshold_to_keep = 0.88;

%% processing loop

if ~exist('ALLEEG','var'); eeglab; end
pop_editoptions( 'option_storedisk', 0, 'option_savetwofiles', 1, 'option_saveversion6', 0, 'option_single', 0, 'option_memmapdata', 0, 'option_eegobject', 0, 'option_computeica', 1, 'option_scaleicarms', 1, 'option_rememberfolder', 1, 'option_donotusetoolboxes', 0, 'option_checkversion', 1, 'option_chat', 1);


for subject = subjects
   
   disp(['Subject #' num2str(subject)]);
   
   % Right now (15-02-2019) there's a bug in mobilab export
   % which shifts marker. Load_xdf is unaffected which is why
   % we use this exclusively for now.
   input_filepath = [process_config.study_folder process_config.raw_data_folder num2str(subject)];
   output_filepath = [process_config.study_folder process_config.raw_EEGLAB_data_folder num2str(subject)];
   
   convert_xdf_to_set(ALLEEG,[],input_filepath,output_filepath);
   
   % This merges all EEG data files into one big file
   input_filepath = [process_config.study_folder process_config.raw_EEGLAB_data_folder num2str(subject)];
   output_filepath = [process_config.study_folder process_config.raw_EEGLAB_data_folder num2str(subject)];
   
   STUDY = []; CURRENTSTUDY = 0; ALLEEG = []; EEG=[]; CURRENTSET=[];
   
   EEG = pop_loadset('filename', strcat(process_config.filenames, '_EEG.set'), 'filepath', input_filepath);
   [ALLEEG, EEG, CURRENTSET] = pop_newset(ALLEEG, EEG, 0,'study',0);
   
   % merges all files currently loaded in EEGLAB into one file and stores
   % the original filenames in EEG.etc.appended_files
   [ALLEEG, EEG_merged, CURRENTSET] = bemobil_merge(ALLEEG,EEG,CURRENTSET,1:length(ALLEEG), process_config.merged_filename, output_filepath);
   
   [ALLEEG, EEG_AMICA_final, CURRENTSET] = bemobil_process_all_AMICA(ALLEEG, EEG_merged, CURRENTSET, subject, process_config);
   
end
