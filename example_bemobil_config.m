%% ONLY CHANGE THESE PARTS!

clear bemobil_config

%% General Setup
bemobil_config.study_folder = 'path_to_study_folder\data\';

bemobil_config.filename_prefix = 'sub_';

bemobil_config.filenames = {'condition1' 'condition2'};
% custom scripts can be applied to the datasets right before splitting the complete MOBI set up. e.g. for parsing
% events, removing breaks and pre/post experiment segments. These are your own scripts! NEED to be taking in three
% parameters: EEGin, filename, filepath and give back one: EEGout 

% There needs to be one entry per filename! leave empty entries if you dont have scripts to apply here
bemobil_config.MOBI_functions = {'parse_events_remove_breaks' 'parse_events_remove_breaks'};
bemobil_config.MOBI_functions = {'' ''}; % if no such scripts exist

% which data streams are in your data? be as specific as possible (the code uses regexp)

% these data streams should be loaded from the xdf file but no mocap processing happens (like EEG and other types)
bemobil_config.unprocessed_data_streams = {'eegstream','eyetracking'}; 
% events of these streams should be exported
bemobil_config.event_streams = {'markers'}; 
% these streams should be processed as rigid body streams containing 3 dof position and 3 dof orientation data (e.g. derivatives and filters applied)
bemobil_config.rigidbody_streams = {'mocap_head','mocap_hand'}; 

% foldernames (NEED to have a filesep at the end, sorry!) and filenames
bemobil_config.raw_data_folder = '0_raw-data\';
bemobil_config.mobilab_data_folder = '1_mobilab-data\';
bemobil_config.raw_EEGLAB_data_folder = '2_basic-EEGLAB\';
bemobil_config.spatial_filters_folder = '3_spatial-filters\';
bemobil_config.spatial_filters_folder_AMICA = '3-1_AMICA\';
bemobil_config.spatial_filters_folder_SSD = '3-2_SSD\';
bemobil_config.single_subject_analysis_folder = '4_single-subject-analysis\';

bemobil_config.merged_filename = 'merged.set';
bemobil_config.preprocessed_filename = 'preprocessed.set';
bemobil_config.interpolated_avRef_filename = 'interpolated_avRef.set';
bemobil_config.filtered_filename = 'filtered.set';
bemobil_config.amica_raw_filename_output = 'postAMICA_raw.set';
bemobil_config.amica_filename_output = 'postAMICA_cleaned.set';
bemobil_config.warped_dipfitted_filename = 'warped_dipfitted.set';
bemobil_config.copy_weights_interpolate_avRef_filename = 'interp_avRef_ICA.set';
bemobil_config.single_subject_cleaned_ICA_filename = 'cleaned_with_ICA.set';
bemobil_config.ssd_frontal_parietal_filename = 'ssd_frontal_parietal.set';

%% Processing Setup 

% enter channels that you did not use at all (e.g. with our custom MoBI 160 chan layout, only 157 chans are used), leave
% empty, if all channels are used
% process_config.channels_to_remove = {'N29' 'N30' 'N31'};
bemobil_config.channels_to_remove = [];

% enter EOG channel names here:
% bemobil_config.eog_channels  = {'G16' 'G32'};
bemobil_config.eog_channels  = {'VEOG'};

% if you add a channel here it needs to have a location as well. this means a new channel will be created and the old reference will be back in the dataset
bemobil_config.ref_channel  = 'Z7(REF)'; 

bemobil_config.rename_channels = []; % it's possible to rename single channels here if needed. Enter matrices of channel names (from->to)

% resample frequency during preprocessing (leave empty if you resample before, or your data is already correctly
% sampled)
bemobil_config.resample_freq = 250; 
% bemobil_config.resample_freq = []; 

% automatic channel cleaning:
% this is the minimal correlation a channel has to have with its own reconstruction based on interpolation, see
% clean_artifacts. 0.8 seems to be reasonable
bemobil_config.chancorr_crit = 0.8;

% channel locations: leave this empty if you have standard channel names that should use standard 10-20 locations,
% otherwise every dataset needs to have a channel locations file in the raw_data folder, and the chanloc file needs to
% have the correct participant prefix!

% bemobil_config.channel_locations_filename = 'channel_locations.elc';
% bemobil_config.channel_locations_filename = [];

% for warping the electrode locations to the standard 10-20 locations (leave
% empty if using standard locations)
% bemobil_config.warping_channel_names = {3,'FTT9h';45,'FTT10h';84,'AFz';87,'Cz'};
bemobil_config.warping_channel_names = [];

% ZapLine to reduce line noise frequencies. You can enter more than one frequency if you have more noise (like from
% lights or TVs) and know the frequency. Leave empty if no noise is present (haha).
bemobil_config.linefreqs = [50]; 
% use adaptive detector to determine the amount of removal by entering empty here (recommended), otherwise set a
% predefined removal, e.g. removing 3 components (this does NOT reduce data rank!)
bemobil_config.zapline_n_remove = [];
bemobil_config.zapline_plot = 1;

% filter for AMICA:
% See Klug & Gramann (2020) for an investigation of filter effect on AMICA -> 1.25 Hz should be a good compromise if you
% don't know how much movement exists, otherwise even higher may be good, up to 2Hz, and you need to subtract 0.25 to
% obtain the correct cutoff value for a filter order of 1650
bemobil_config.filter_lowCutoffFreqAMICA = 1.25; % 1.25 is 1Hz cutoff!
bemobil_config.filter_AMICA_highPassOrder = 1650; % was used by Klug & Gramann (2020)
bemobil_config.AMICA_autoreject = 1; % uses automatic rejection method of AMICA. no time-cleaning (manual or automatic) is needed then!

%% Special Processing Parameters
% everything from here is according to the general pipeline, changes not recommended 

% mocap processing

bemobil_config.mocap_lowpass = 6;
bemobil_config.rigidbody_derivatives = 2;

%%% AMICA

% on some PCs AMICA may crash before the first iteration if the number of threads and the amount the data does not suit
% the algorithm. Jason Palmer has been informed, but no fix so far. just roll with it. if you see the first iteration
% working there won't be any further crashes. in this case just press "close program" or the like and the
% bemobil_spatial_filter algorithm will AUTOMATICALLY reduce the number of threads and start AMICA again. this way you
% will always have the maximum number of threads that should be used for AMICA. check in the task manager how many
% threads you have theoretically available and think how much computing power you want to devote for AMICA. 

% 4 threads are most effective for single subject speed, more threads don't really shorten the calculation time much.
% best efficiency is using just 1 thread and have as many matlab instances open as possible (limited by the CPU usage).
% Remember your RAM limit in this case.

bemobil_config.filter_highCutoffFreqAMICA = [];
bemobil_config.filter_AMICA_lowPassOrder = [];
bemobil_config.max_threads = 4;
bemobil_config.num_models = 1;

% dipfit settings
bemobil_config.residualVariance_threshold = 100;
bemobil_config.do_remove_outside_head = 'off';
bemobil_config.number_of_dipoles = 1;

% IC_label settings
% 'default' classifier did not lead to good classification of muscles (see Klug & Gramann (2020)), 'lite' was better
% overall.
bemobil_config.iclabel_classifier = 'lite';

% 'Brain', 'Muscle', 'Eye', 'Heart', 'Line Noise', 'Channel Noise', 'Other'
bemobil_config.iclabel_classes = [1];

% if the threshold is set to -1, the popularity classifier is used (i.e. every IC gets the class with the highest
% probability), if it is set to a value, the summed score of the iclabel_classes must be higher than this threshold to
% keep an IC. Must be in the [0 1] range!
bemobil_config.iclabel_threshold = -1; 


% for SSD analysis (ignore if you don't know)
bemobil_config.frontal_channames = {'Fz','FCz','F1','F2'};
bemobil_config.parietal_channames = {'Pz','P1','P2','P3','P4'};

% SSD analysis
bemobil_config.ssd_freq_theta = [4 7; % signal band
	1 10; % noise bandbass (outer edges)
	2.5 8.5]; % noise bandstop (inner edges)
bemobil_config.ssd_freq_alpha = [8 13; % signal band
	5 16; % noise bandbass (outer edges)
	6.5 14.5]; % noise bandstop (inner edges)

