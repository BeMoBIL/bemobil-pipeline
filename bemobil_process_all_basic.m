function [ALLEEG, EEG_interp_avRef, CURRENTSET] = bemobil_process_all_basic(subject, bemobil_config, ALLEEG, EEG_merged, CURRENTSET, force_recompute)

if ~exist('force_recompute','var') || isempty(force_recompute)
    force_recompute = 0;
end


disp(['Subject #' num2str(subject)]);

% get rid of memory mapped object storage and make sure double spacing and matlab save version 7 is used (for files
% larger than 2gb)
% mobilab uses memory mapped files which is why this needs to be set several times throughout the processing
try
    pop_editoptions( 'option_saveversion6', 0, 'option_single', 0, 'option_memmapdata', 0);
catch
    warning('Could NOT edit EEGLAB memory options!!');
end


filepath = fullfile(bemobil_config.study_folder, bemobil_config.raw_EEGLAB_data_folder, [bemobil_config.filename_prefix num2str(subject)]);


% check if the whole script has been running already
if ~force_recompute
    try
        EEG = pop_loadset('filename', [bemobil_config.filename_prefix num2str(subject) '_'...
            bemobil_config.interpolated_avRef_filename], 'filepath', filepath);
        [ALLEEG, EEG_interp_avRef, CURRENTSET] = pop_newset(ALLEEG, EEG, 0,'study',0);
        
        warning('Old interpolated file already existed, using that file!')
        
        return
    catch
        disp('...failed.')
    end
end


% preprocess
if ~force_recompute
    try
        
        EEG_preprocessed = pop_loadset('filename', [bemobil_config.filename_prefix num2str(subject) '_'...
            bemobil_config.preprocessed_filename], 'filepath', filepath);
        
        warning('Old preprocessed file already existed, using that file!')
        
        % save RAM
        clear EEG_merged
    catch
        disp('...failed.')
    end
end

if ~exist('EEG_preprocessed','var')
    
    if ~isempty(bemobil_config.channel_locations_filename)
        channel_locations_filepath = fullfile(bemobil_config.study_folder, bemobil_config.raw_data_folder,...
            [bemobil_config.filename_prefix num2str(subject)], [bemobil_config.filename_prefix num2str(subject) '_'...
            bemobil_config.channel_locations_filename]);
    else
        channel_locations_filepath = [];
    end
    
    % preprocessing: enter chanlocs, remove unused channels, declare EOG, resample
    [ALLEEG, EEG_preprocessed, CURRENTSET] = bemobil_preprocess(ALLEEG, EEG_merged, CURRENTSET, channel_locations_filepath,...
        bemobil_config.channels_to_remove, bemobil_config.eog_channels, bemobil_config.resample_freq,...
        [bemobil_config.filename_prefix num2str(subject) '_' bemobil_config.preprocessed_filename], filepath,...
        bemobil_config.rename_channels, bemobil_config.ref_channel, bemobil_config.zaplineConfig);
    
    disp('Preprocessing done!')
    
    % save RAM
    clear EEG_merged
end

%% bad channel removal with clean_artifacts

EEG = EEG_preprocessed;

%% detect bad channels

[chans_to_interp, plothandle] = bemobil_detect_bad_channels(EEG, ALLEEG, CURRENTSET, bemobil_config.chancorr_crit,...
    bemobil_config.chan_max_broken_time);

%% save fig of bad channels

savefig(plothandle,fullfile(filepath,[bemobil_config.filename_prefix num2str(subject) '_preprocessed_bad_channels.fig']))
print(plothandle,fullfile(filepath,[bemobil_config.filename_prefix num2str(subject) '_preprocessed_bad_channels.png']),'-dpng')
close

%% do the actual interpolation and average referencing (reference is not considering EOGs)

disp('Interpolating bad channels and compute final average reference, ignoring EOG channels...')
[ALLEEG, EEG_interp_avRef, CURRENTSET] = bemobil_interp_avref( EEG_preprocessed , ALLEEG, CURRENTSET, chans_to_interp,...
    [bemobil_config.filename_prefix num2str(subject) '_' bemobil_config.interpolated_avRef_filename], filepath);

%% plot interpolated filtered, for analytics

disp('Filtering data only for plotting!')
EEG = pop_eegfiltnew(EEG_interp_avRef, 'locutoff',0.5);

%%
bemobil_plot_data_chunks(EEG, filepath, [bemobil_config.filename_prefix num2str(subject) '_preprocessed_interpolated_channels'])

disp('All basic EEG processing done.')