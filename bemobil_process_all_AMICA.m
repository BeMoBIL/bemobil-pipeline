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

% get rid of memory mapped object storage and make sure double spacing and matlab save version 7 is used (for files
% larger than 2gb)
% mobilab uses memory mapped files which is why this needs to be set several times throughout the processing
try
    pop_editoptions( 'option_saveversion6', 0, 'option_single', 0, 'option_memmapdata', 0);
catch
    warning('Could NOT edit EEGLAB memory options!!');
end

if ~exist('force_recompute','var')
    force_recompute = false;
end
if force_recompute
    warning('RECOMPUTING OLD FILES IF FOUND!!!')
end

% check if the entire processing was done already
output_filepath = fullfile(bemobil_config.study_folder, bemobil_config.single_subject_analysis_folder, [bemobil_config.filename_prefix num2str(subject)]);

if ~force_recompute
    try
        
        EEG = pop_loadset('filename', [bemobil_config.filename_prefix num2str(subject) '_'...
            bemobil_config.copy_weights_interpolate_avRef_filename], 'filepath', output_filepath);
        [ALLEEG, EEG_single_subject_final, CURRENTSET] = pop_newset(ALLEEG, EEG, 0,'study',0);
        
        warning('Old single subject file with interpolated channels, avref, and AMICA data already existed, using that file!')
        
    end
end

if ~exist('EEG_single_subject_final','var')
    
    output_filepath = fullfile(bemobil_config.study_folder, bemobil_config.spatial_filters_folder,...
        bemobil_config.spatial_filters_folder_AMICA, [bemobil_config.filename_prefix num2str(subject)]);
    
    % check if the part was done already
    if ~force_recompute
        try
            
            EEG = pop_loadset('filename', [bemobil_config.filename_prefix num2str(subject) '_'...
                bemobil_config.amica_filename_output], 'filepath', output_filepath);
            [ALLEEG, EEG_AMICA_cleaned, CURRENTSET] = pop_newset(ALLEEG, EEG, 0,'study',0);
            
            warning('Old AMICA file already existed, using that file!')
            
        end
    end
    
    if ~exist('EEG_AMICA_cleaned','var')
        
        %% highpass filter for AMICA
        
        [ALLEEG, EEG_filtered_for_AMICA, CURRENTSET] = bemobil_filter(ALLEEG, EEG_interp_avRef, CURRENTSET,...
            bemobil_config.filter_lowCutoffFreqAMICA, bemobil_config.filter_highCutoffFreqAMICA,...
            [], [], bemobil_config.filter_AMICA_highPassOrder, bemobil_config.filter_AMICA_lowPassOrder);
        
        % save RAM and disk space of ICA results, since events are irrelevant here and in mobi datasets can be a lot
        EEG_filtered_for_AMICA.event = [];
        EEG_filtered_for_AMICA.urevent = [];
        
        %% AMICA
        % running signal decomposition with automatic rejection is recommended
        
        %         turns out the data rank reduction due to bridges does not really solve the issue of a high D value when
        %         computing AMICA, and since PCA before ICA was shown the be problematic, I removed it again: Artoni, F.,
        %         Delorme, A., & Makeig, S. (2018) Applying dimension reduction to EEG data by Principal Component Analysis
        %         reduces the quality of its subsequent Independent Component decomposition. Neuroimage, 175, 176–187.
        
        %         data_rank = EEG_filtered_for_AMICA.etc.rank;
        %         [rank_reduction_of_bridges,EEG_filtered_for_AMICA] = bemobil_find_gel_bridges(EEG_filtered_for_AMICA,0.98);
        %         data_rank = data_rank - rank_reduction_of_bridges;
        
        [ALLEEG, EEG_AMICA_cleaned, CURRENTSET] = bemobil_signal_decomposition(ALLEEG, EEG_filtered_for_AMICA, ...
            CURRENTSET, true, bemobil_config.num_models, bemobil_config.max_threads, EEG_filtered_for_AMICA.etc.rank, [], ...
            [bemobil_config.filename_prefix num2str(subject) '_' bemobil_config.amica_filename_output], output_filepath,...
            bemobil_config.AMICA_autoreject, bemobil_config.AMICA_n_rej, bemobil_config.AMICA_reject_sigma_threshold);
        
        % add information about AMICA autorejected time points
        
        sample_mask = EEG_AMICA_cleaned.etc.spatial_filter.AMICAmods.Lt == 0;
        EEG_AMICA_cleaned.etc.bad_samples = sample_mask;
        EEG_AMICA_cleaned.etc.bad_samples_percent = sum(EEG_AMICA_cleaned.etc.bad_samples) / length(EEG_AMICA_cleaned.etc.bad_samples) * 100;
        
        % find latency of regions
        remove_data_intervals = reshape(find(diff([false sample_mask false])),2,[])';
        remove_data_intervals(:,2) = remove_data_intervals(:,2)-1;
        EEG_AMICA_cleaned.etc.remove_data_intervals = remove_data_intervals;
        
        % save again
        EEG_AMICA_cleaned = pop_saveset( EEG_AMICA_cleaned,...
            'filename',[bemobil_config.filename_prefix num2str(subject) '_' bemobil_config.amica_filename_output],...
            'filepath', output_filepath);
        disp('...done');
        
        bemobil_plot_amica_autoreject(EEG_AMICA_cleaned, output_filepath, [bemobil_config.filename_prefix num2str(subject) '_AMICA_autoreject'])
        
        % save RAM
        clear EEG_filtered_for_AMICA
        
    end
    
    
    %% Warping of locations and dipole fitting, plus runing ICLabel
    % renames the specified channels, warps the chanlocs on a standard head model and fits dipoles for
    % each IC below the threshold of residual variance
    
    output_filepath = fullfile(bemobil_config.study_folder, bemobil_config.spatial_filters_folder,...
        bemobil_config.spatial_filters_folder_AMICA, [bemobil_config.filename_prefix num2str(subject)]);
    
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
    
    % plot all AMICA results
    pop_topoplot(EEG, 0, [1:size(EEG.icaweights,1)] ,EEG.filename,[] ,0,'electrodes','off'); allICfighandle = gcf;
    print(allICfighandle,fullfile(output_filepath,[bemobil_config.filename_prefix num2str(subject) '_all_ICs.png']),'-dpng')
    close
    
    %% Copy the spatial filter data into the raw full data set for further single subject processing
    
    output_filepath = fullfile(bemobil_config.study_folder, bemobil_config.single_subject_analysis_folder,...
        [bemobil_config.filename_prefix num2str(subject)]);
    
    disp('Copying all information into full length dataset for single subject processing...');
    [ALLEEG, EEG_single_subject_copied, CURRENTSET] = bemobil_copy_spatial_filter(EEG_interp_avRef, ALLEEG, CURRENTSET,...
        EEG_AMICA_final, [bemobil_config.filename_prefix num2str(subject) '_'...
        bemobil_config.copy_weights_interpolate_avRef_filename], output_filepath);
    
    % save RAM
    clear EEG_interp_avRef EEG_AMICA_final
    %% clean with IClabel
    
    disp('Cleaning data with ICLabel')
    
    % clean now, save files and figs
    [ALLEEG, EEG_single_subject_final, CURRENTSET, ICs_keep, ICs_throw] = bemobil_clean_with_iclabel( EEG_single_subject_copied ,...
        ALLEEG, CURRENTSET, bemobil_config.iclabel_classifier,...
        bemobil_config.iclabel_classes, bemobil_config.iclabel_threshold,...
        [ bemobil_config.filename_prefix num2str(subject) '_' bemobil_config.single_subject_cleaned_ICA_filename],output_filepath);
    
    disp('...done.')
    
    %% plot interpolated filtered, for analytics
    
    disp('Filtering data only for plotting!')
    EEG = pop_eegfiltnew(EEG_single_subject_final, 'locutoff',0.5);
    
    bemobil_plot_data_chunks(EEG, output_filepath, [bemobil_config.filename_prefix num2str(subject) '_cleaned'])

    
end

disp('Entire AMICA processing done!');
