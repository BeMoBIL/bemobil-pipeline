% bemobil_check_config - Checks if all fields of the config are present. If fields are missing, they will be shown, and
% in case a default for this field exists it will be entered. If no default exists, the filed is considered mandatory
% and an error is thrown.
%
% Usage:
%   >>  [bemobil_config] = bemobil_check_config(bemobil_config)
% 
% Inputs:
%   bemobil_config              - current bemobil_config struct
%
% Outputs:
%   bemobil_config              - checked bemobil_config struct with default filled if missing
%
% See also:
%   example_bemobil_config
%
% Authors: Marius Klug, 2021

function bemobil_config = bemobil_check_config(bemobil_config)

if nargin == 0
    help bemobil_check_config
    return
end

% remove excess fieldtrip paths in case they are there
ftpath = fileparts(which('ft_defaults'));
disp('Removing fieldtrip paths that would otherwise replace MATLAB toolboxes...')
rmpath(fullfile(ftpath, 'external', 'signal'))
rmpath(fullfile(ftpath, 'external', 'stats'))
rmpath(fullfile(ftpath, 'external', 'image'))

all_config_fields = {
    'study_folder', {}
    'filename_prefix', {'sub-'}
    'source_data_folder', {['0_source-data' filesep]}
    'bids_data_folder', {['1_BIDS-data' filesep]}
    'raw_EEGLAB_data_folder', {['2_raw-EEGLAB' filesep]}
    'EEG_preprocessing_data_folder', {['3_EEG-preprocessing' filesep]}
    'spatial_filters_folder', {['4_spatial-filters' filesep]}
    'spatial_filters_folder_AMICA', {['4-1_AMICA' filesep]}
    'single_subject_analysis_folder', {['5_single-subject-EEG-analysis' filesep]}
    'motion_analysis_folder', {['6_single-subject-motion-analysis' filesep]}
    'merged_filename', {'merged_EEG.set'}
    'basic_prepared_filename', {'basic_prepared.set'}
    'preprocessed_filename', {'preprocessed.set'}
    'filtered_filename', {'filtered.set'}
    'amica_filename_output', {'AMICA.set'}
    'dipfitted_filename', {'dipfitted.set'}
    'preprocessed_and_ICA_filename', {'preprocessed_and_ICA.set'}
    'single_subject_cleaned_ICA_filename', {'cleaned_with_ICA.set'}
    'merged_motion_filename', {'merged_MOTION.set'}
    'processed_motion_filename', {'motion_processed.set'}
    'channels_to_remove',{[]}
    'eog_channels',{[]}
    'ref_channel',{[]}
    'rename_channels',{[]}
    'resample_freq',{250}
    'chancorr_crit',{0.8}
    'chan_max_broken_time',{0.5}
    'chan_detect_num_iter',{20}
    'chan_detected_fraction_threshold',{0.5}
    'flatline_crit',{'off'}
    'line_noise_crit',{'off'}
    'num_chan_rej_max_target',{1/5}
    'channel_locations_filename',{[]}
    'zaplineConfig',{struct('noisefreqs',[])}
    'filter_lowCutoffFreqAMICA',{1.75}
    'filter_AMICA_highPassOrder',{1650}
    'filter_highCutoffFreqAMICA',{[]}
    'filter_AMICA_lowPassOrder',{[]}
    'num_models',{1}
    'AMICA_autoreject',{1}
    'AMICA_n_rej',{10}
    'AMICA_reject_sigma_threshold',{3}
    'AMICA_max_iter',{2000}
    'max_threads',{4}
    'warping_channel_names',{[]}
    'residualVariance_threshold',{100}
    'do_remove_outside_head',{'off'}
    'number_of_dipoles',{1}
    'iclabel_classifier',{'default'}
    'iclabel_classes',{[1 2 4 5 6 7]}
    'iclabel_threshold',{-1}
    'final_filter_lower_edge',{0.2}
    'final_filter_higher_edge',{[]}
    'lowpass_motion',{6}
    'lowpass_motion_after_derivative',{18}
    };

%% check fields

disp('Checking config fields...')

for i_field = 1:length(all_config_fields)
    
    if ~isfield(bemobil_config,all_config_fields{i_field,1})
        warning(['Field not present in bemobil_config: ''' all_config_fields{i_field,1} ''''])
        
        if ~isempty(all_config_fields{i_field,2})
            this_default = all_config_fields{i_field,2};
            disp('Using default:')
            disp(this_default{1})
            bemobil_config.(all_config_fields{i_field,1}) = this_default{1};
        else
            error(['Mandatory field missing in bemobil_config: bemobil_config.' all_config_fields{i_field,1}])
        end
    end
end

disp('All fields checked!')

try 
    pop_editoptions('option_saveversion6', 0, 'option_single', 0, 'option_memmapdata', 0, 'option_savetwofiles', 1, 'option_storedisk', 0);
catch
    warning('Could NOT edit EEGLAB memory options!!'); 
end
