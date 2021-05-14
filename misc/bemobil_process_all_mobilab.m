% bemobil_process_all_mobilab - wrapper function that incorporates all necessary processing steps from raw .xdf to .set
% files in EEGLAB. Data is being loaded into mobilab, rigidbody motion streams are processed (filtered, transformed into
% euler angles, derived) and data and marker streams are then exported to EEGLAB containing EEG and all other kinds of
% data. The dataset has the suffix '_MoBI'. This dataset is then split into individual channel types (e.g. 'EEG',
% 'MOTION', 'EYE', 'OTHER'), and subsequently all EEG files (from several raw .xdf files) will be merged into one large
% EEG file for this participant, which can then be used for further processing (e.g. with bemobil_process_all_AMICA)
%
% The intermediate files are stored on the disk.
%
% Usage:
%   >>  [ALLEEG, EEG_merged, CURRENTSET] = bemobil_process_all_mobilab(subject, bemobil_config, ALLEEG, CURRENTSET, mobilab)
%
% Inputs:
%   subject                   - subject number of the current subject (necessary for filepaths and storage)
%   bemobil_config            - configuration struct with all necessary information. See EEG_processing_example file
%                                that comes with this function!
%   ALLEEG                    - complete EEGLAB data set structure
%   CURRENTSET                - index of current EEGLAB EEG structure within ALLEEG
%	mobilab					  - container for the mobilab application. execute "runmobilab" before this script to get it
%
% Outputs:
%   ALLEEG                    - complete EEGLAB data set structure
%   EEG_interp_avRef		  - preprocessed and interpolated EEGLAB EEG structure, contains EEG datasets of all conditions.
%   CURRENTSET                - index of current EEGLAB EEG structure within ALLEEG
%
%   .set data file of EEGLAB EEG structures are stored on disk according to their names in the bemobil_config
%
% See also:
%   EEGLAB
%
% Authors: Marius Klug, 2019

function [ALLEEG, EEG_interp_avRef, CURRENTSET] = bemobil_process_all_mobilab(subject, bemobil_config, ALLEEG, CURRENTSET, mobilab, force_recompute)

if ~exist('force_recompute','var')
    force_recompute = false;
end

if force_recompute
    warning('RECOMPUTING OLD FILES IF FOUND!!!')
end

disp(['Subject #' num2str(subject)]);

input_filepath = fullfile(bemobil_config.study_folder, bemobil_config.raw_data_folder, [bemobil_config.filename_prefix num2str(subject)]);
output_filepath_mobi = fullfile(bemobil_config.study_folder, bemobil_config.mobilab_data_folder, [bemobil_config.filename_prefix num2str(subject)]);
output_filepath = fullfile(bemobil_config.study_folder, bemobil_config.raw_EEGLAB_data_folder, [bemobil_config.filename_prefix num2str(subject)]);

% get rid of memory mapped object storage and make sure double spacing and matlab save version 7 is used (for files
% larger than 2gb)
% mobilab uses memory mapped files which is why this needs to be set several times throughout the processing
try
    pop_editoptions( 'option_saveversion6', 0, 'option_single', 0, 'option_memmapdata', 0);
catch
    warning('Could NOT edit EEGLAB memory options!!');
end

% make sure EEGLAB has no files other than the ones to be merged
% STUDY = []; CURRENTSTUDY = 0; ALLEEG = []; EEG=[]; CURRENTSET=[];


% check if the whole script has been running already
if ~force_recompute
    try
        
        EEG = pop_loadset('filename', [bemobil_config.filename_prefix num2str(subject) '_'...
            bemobil_config.interpolated_avRef_filename], 'filepath', output_filepath);
        [ALLEEG, EEG_interp_avRef, CURRENTSET] = pop_newset(ALLEEG, EEG, 0,'study',0);
        
        warning('Old interpolated file already existed, using that file!')
        
        return
    end
    
    % check if the full merged EEG file exists
    try
        
        EEG = pop_loadset('filename', [bemobil_config.filename_prefix num2str(subject) '_'...
            bemobil_config.merged_filename], 'filepath', output_filepath);
        [ALLEEG, EEG_merged, CURRENTSET] = pop_newset(ALLEEG, EEG, 0,'study',0);
        
        warning('Old merged file already existed, using that file!')
        
    end
end
%% run MoBILAB

if ~exist('EEG_merged','var') || force_recompute
    
    
    % check if at least the mobilab part was done already and all MoBI
    % files exist
    mobi_files_exist = zeros(size(bemobil_config.filenames));
    for i_filename = 1:length(bemobil_config.filenames)
        
        full_filename = [bemobil_config.filename_prefix num2str(subject) '_' bemobil_config.filenames{i_filename}];
        
        try
            MoBI_EEG = pop_loadset('filename',[full_filename '_MoBI.set'],'filepath', output_filepath);
            mobi_files_exist(i_filename) = 1;
        catch
            mobi_files_exist(i_filename) = 0;
        end
    end
    
    if any(mobi_files_exist) && ~force_recompute
        warning('Old MoBI files already existed, using these files!')
        datasets_to_compute = find(~mobi_files_exist);
        % skipping the mobilab part for these sets
    else
        datasets_to_compute = 1:length(bemobil_config.filenames);
    end
    % doing the mobilab part
    for i_filename = datasets_to_compute
        
        full_filename = [bemobil_config.filename_prefix num2str(subject) '_' bemobil_config.filenames{i_filename}];
        
        %% import from mobilab
        disp(['Importing file: "' full_filename '.xdf" ...']);
        
        input_filepath_mobilab = fullfile(input_filepath, [full_filename '.xdf']);
        
        % this suffix needs to be changed if you redo it and want to keep the old files
        output_filepath_mobilab = fullfile(output_filepath_mobi, [full_filename '_MoBI']);
        
        files_in_mobifolder = dir(output_filepath_mobilab);
        
        if ~isempty(files_in_mobifolder) && force_recompute
            
            % if old files exist, mobilab attempts to zip them, which takes ages, so they are deleted here
            warning('MoBI folder already existed, deleting!')
            
            status = rmdir(output_filepath_mobilab,'s');
            
            if ~status
                % zipping of small text files is fast
                warning('Could not delete all old MoBI files, a zip file with the old logfile and notes will remain.')
            end
        elseif length({files_in_mobifolder.name})>4 && ~force_recompute
            warning('MoBI folder already existed, using old file!')
            continue
        end
        
        % now load
        mobilab.allStreams = dataSourceXDF(input_filepath_mobilab,output_filepath_mobilab);
        all_mobilab_streams = [mobilab.allStreams().item];
        mobilab.refresh
        
        disp('...done.');
        
        
        %% do rigidbody processing
        
        % variable for storing the stream names for exporting later
        processed_rb_stream_names = {};
        
        for i_rigidbody_stream = 1:length(bemobil_config.rigidbody_streams)
            
            % find out the index of this rb stream
            for i_stream = 1:length(all_mobilab_streams)
                
                if regexp(lower(all_mobilab_streams{i_stream}.name),lower(bemobil_config.rigidbody_streams{i_rigidbody_stream}))
                    
                    mobilab_rb_index = i_stream;
                    
                    % process motion channels
                    
                    % quaternion values can flip their signs, which does indicate the same orientation. we need them to be
                    % a smooth curve for filtering, though
                    unflipped = mobilab.allStreams().item{mobilab_rb_index}.unflipSigns();
                    
                    % lowpass filter with the respective frequency (e.g. 6 Hz), zero-lag FIR filter
                    filtered = unflipped.lowpass(bemobil_config.motion_lowpass);
                    
                    % transform quaternion values into Euler angles
                    eulers = filtered.quaternionsToEuler();
                    
                    % compute derivatives (e.g. 2 -> velocity and acceleration)
                    eulers.timeDerivative(bemobil_config.rigidbody_derivatives);
                    
                    % update GUI
                    mobilab.refresh
                    
                    
                    % store names for exporting
                    processed_rb_stream_names{end+1} = eulers.name;
                    
                    % 		child names don't need to be stored since they contain the eulers names and regexp is used later
                    % 		for i_child = 1:bemobil_config.rigidbody_derivatives
                    %
                    % 			processed_rb_stream_names{end+1} = eulers.children{i_child}.name;
                    %
                    % 		end
                    
                    break
                end
                
            end
            
        end
        
        %% find out the index of all data and event streams
        
        all_data_stream_names = [bemobil_config.unprocessed_data_streams,processed_rb_stream_names];
        all_data_stream_indices = [];
        all_event_stream_indices = [];
        
        all_mobilab_streams = [mobilab.allStreams().item];
        
        for i_data_stream = 1:length(all_data_stream_names)
            for i_stream = 1:length(all_mobilab_streams)
                if regexp(lower(all_mobilab_streams{i_stream}.name),lower(all_data_stream_names{i_data_stream}))
                    
                    all_data_stream_indices(end+1) = i_stream;
                    
                end
            end
        end
        
        for i_event_stream = 1:length(bemobil_config.event_streams)
            for i_stream = 1:length(all_mobilab_streams)
                if regexp(lower(all_mobilab_streams{i_stream}.name),lower(bemobil_config.event_streams{i_event_stream}))
                    
                    all_event_stream_indices(end+1) = i_stream;
                    
                end
            end
        end
        
        %% export MoBI dataset
        
        disp('Exporting MoBI dataset. This may take a while!')
        exported_EEG = mobilab.allStreams().export2eeglab(all_data_stream_indices,all_event_stream_indices);
        
        % clear RAM
        mobilab.allStreams = [];
        
        mkdir(output_filepath)
        
        % get rid of memory mapped object storage and make sure double spacing and matlab save version 7 is used (for files
        % larger than 2gb)
        % mobilab uses memory mapped files which is why this needs to be set several times throughout the processing
        try
            pop_editoptions( 'option_saveversion6', 0, 'option_single', 0, 'option_memmapdata', 0);
        catch
            warning('Could NOT edit EEGLAB memory options!!');
        end
        
        pop_saveset( exported_EEG, 'filename',[full_filename '_MoBI'],'filepath', output_filepath);
        disp('...done');
        clear exported_EEG
        
    end % mobilab loops
    
    %% load and merge MoBI files
    % get rid of memory mapped object storage and make sure double spacing and matlab save version 7 is used (for files
    % larger than 2gb)
    % mobilab uses memory mapped files which is why this needs to be set several times throughout the processing
    try
        pop_editoptions( 'option_saveversion6', 0, 'option_single', 0, 'option_memmapdata', 0);
    catch
        warning('Could NOT edit EEGLAB memory options!!');
    end
    
    for i_filename = 1:length(bemobil_config.filenames)
        
        full_filename = [bemobil_config.filename_prefix num2str(subject) '_' bemobil_config.filenames{i_filename}];
        
        MoBI_EEG = pop_loadset('filename',[full_filename '_MoBI.set'],'filepath', output_filepath);
        
        if isfield(bemobil_config, 'MOBI_functions')
            if ~isempty(bemobil_config.MOBI_functions{i_filename})
                % this allows for custom functions to happen before splitting the MOBI dataset
                MoBI_EEG = feval(bemobil_config.MOBI_functions{i_filename}, MoBI_EEG, full_filename, output_filepath);
            end
            MoBI_EEG.etc.applied_MoBI_script = bemobil_config.MOBI_functions{i_filename};
        end
        
        % split MoBI dataset into unique channel types
        [~, ~, ~, EEG_split_sets] = bemobil_split_MoBI_set(ALLEEG, MoBI_EEG, CURRENTSET);
        
        % clear RAM
        STUDY = []; CURRENTSTUDY = 0; ALLEEG = [];  CURRENTSET=[];
        MoBI_EEG=[];
        
        for i_splitset = 1:length(EEG_split_sets)
            pop_saveset( EEG_split_sets(i_splitset), 'filename',[full_filename '_'...
                EEG_split_sets(i_splitset).chanlocs(1).type],'filepath', output_filepath);
            disp('...done');
        end
        
        % clear RAM
        EEG_split_sets=[];
        
    end
    
    % This merges all EEG data files into one big file
    input_filepath = [bemobil_config.study_folder bemobil_config.raw_EEGLAB_data_folder bemobil_config.filename_prefix num2str(subject)];
    output_filepath = [bemobil_config.study_folder bemobil_config.raw_EEGLAB_data_folder bemobil_config.filename_prefix num2str(subject)];
    
    % make sure EEGLAB has no files other than the ones to be merged
    STUDY = []; CURRENTSTUDY = 0; ALLEEG = []; EEG=[]; CURRENTSET=[];
    
    EEG = pop_loadset('filename', strcat(strcat(bemobil_config.filename_prefix, num2str(subject), '_',...
        bemobil_config.filenames,'_EEG.set')), 'filepath', input_filepath);
    [ALLEEG, EEG, CURRENTSET] = pop_newset(ALLEEG, EEG, 0,'study',0);
    
    % merges all files currently loaded in EEGLAB into one file and stores
    % the original filenames in EEG.etc.appended_files
    % added LG 02.10.2019
    % in case there is only 1 file do not merge
    if size(ALLEEG,2)>1
        [ALLEEG, EEG_merged, CURRENTSET] = bemobil_merge(ALLEEG,EEG,CURRENTSET,1:length(ALLEEG),...
            [bemobil_config.filename_prefix num2str(subject) '_' bemobil_config.merged_filename], output_filepath);
    else
        EEG_merged = EEG;
    end
    
    
    disp('Entire mobilab loading, processing, and exporting done!')
    disp('You can start the EEGLAB processing now, using the merged dataset.')
    
    clear mobilab
    STUDY = []; CURRENTSTUDY = 0; ALLEEG = []; EEG=[]; CURRENTSET=[];
    %%
end % loading and splitting data to get merged EEG 

%% preprocess

% get rid of memory mapped object storage and make sure double spacing and matlab save version 7 is used (for files
% larger than 2gb)
% mobilab uses memory mapped files which is why this needs to be set several times throughout the processing
try
    pop_editoptions( 'option_saveversion6', 0, 'option_single', 0, 'option_memmapdata', 0);
catch
    warning('Could NOT edit EEGLAB memory options!!');
end

try
    
    EEG_preprocessed = pop_loadset('filename', [bemobil_config.filename_prefix num2str(subject) '_'...
        bemobil_config.preprocessed_filename], 'filepath', output_filepath);
    
    warning('Old preprocessed file already existed, using that file!')
    
    % save RAM
    clear EEG_merged
    
end

if ~exist('EEG_preprocessed','var')
    
    output_filepath = fullfile(bemobil_config.study_folder, bemobil_config.raw_EEGLAB_data_folder, [bemobil_config.filename_prefix num2str(subject)]);
    if ~isempty(bemobil_config.channel_locations_filename)
        channel_locations_filepath = fullfile(bemobil_config.study_folder, bemobil_config.raw_data_folder,...
            [bemobil_config.filename_prefix num2str(subject)], [bemobil_config.filename_prefix num2str(subject) '_'...
            bemobil_config.channel_locations_filename]);
    else
        channel_locations_filepath = [];
    end
    
    % preprocessing: enter chanlocs, remove unused channels, declare EOG, resample
    [ALLEEG, EEG_preprocessed, CURRENTSET] = bemobil_process_EEG_basics(ALLEEG, EEG_merged, CURRENTSET, channel_locations_filepath,...
        bemobil_config.channels_to_remove, bemobil_config.eog_channels, bemobil_config.resample_freq,...
        [bemobil_config.filename_prefix num2str(subject) '_' bemobil_config.preprocessed_filename], output_filepath,...
        bemobil_config.rename_channels, bemobil_config.ref_channel, bemobil_config.linefreqs, bemobil_config.zapline_n_remove,...
        bemobil_config.zapline_plot);
    
    disp('Preprocessing done!')
    
    % save RAM
    clear EEG_merged
end
%% bad channel removal with clean_artifacts

EEG = EEG_preprocessed;

%%

% compute average reference before finding bad channels 
[ALLEEG, EEG, CURRENTSET] = bemobil_avref( EEG , ALLEEG, CURRENTSET);

disp('Detecting bad channels...')

% remove bad channels, use default values of clean_artifacts, but specify just in case they may change
[EEG_chan_removed,EEG_highpass,~,chans_to_interp] = clean_artifacts(EEG,...
    'burst_crit','off','window_crit','off',...
    'chancorr_crit',bemobil_config.chancorr_crit,'line_crit',4,'highpass_band',[0.25 0.75],'flatline_crit','off');
chans_to_interp = find(chans_to_interp); % transform logical array to indices

disp('Detected bad channels: ')
disp({EEG.chanlocs(chans_to_interp).labels})

% take EOG out of the channels to interpolate: EOG is very likely to be different from the others but rightfully so
disp('Ignoring EOG channels for interpolation:')

disp({EEG.chanlocs(chans_to_interp(strcmp({EEG.chanlocs(chans_to_interp).type},'EOG'))).labels})
chans_to_interp(strcmp({EEG.chanlocs(chans_to_interp).type},'EOG'))=[];

disp('Final bad channels: ')
disp({EEG.chanlocs(chans_to_interp).labels})

EEG_chan_removed = pop_select( EEG_highpass,'nochannel',chans_to_interp);
EEG_chan_removed.etc.clean_channel_mask = true(EEG_highpass.nbchan,1);
EEG_chan_removed.etc.clean_channel_mask(chans_to_interp) = deal(0);

% display 1/10 of the data in the middle (save disk space when saving figure)
vis_artifacts(EEG_chan_removed,EEG,'show_events',0,'time_subset',...
    [round(EEG.times(end)/2) round(EEG.times(end)/2+round(EEG.times(end)/10))]/1000);

%%
set(gcf, 'Position', get(0,'screensize'))

savefig(gcf,fullfile(output_filepath,[bemobil_config.filename_prefix num2str(subject) '_preprocessed_bad_channels.fig']))
print(gcf,fullfile(output_filepath,[bemobil_config.filename_prefix num2str(subject) '_preprocessed_bad_channels.png']),'-dpng')
close

%% do the interpolation and average referencing (reference is not considering EOGs)

disp('Interpolating bad channels and compute final average reference, ignoring EOG channels...')
[ALLEEG, EEG_interp_avRef, CURRENTSET] = bemobil_interp_avref( EEG_preprocessed , ALLEEG, CURRENTSET, chans_to_interp,...
    [bemobil_config.filename_prefix num2str(subject) '_' bemobil_config.interpolated_avRef_filename], output_filepath);

%% plot interpolated filtered, for analytics

EEG = pop_eegfiltnew(EEG_interp_avRef, 'locutoff',0.5);
EEG.filename = [bemobil_config.filename_prefix num2str(3) '_interpAveRef_highpass'];

vis_artifacts(EEG,EEG,'show_events',0,'time_subset',...
    [round(EEG.times(end)/2) round(EEG.times(end)/2+round(EEG.times(end)/10))]/1000);

%%
set(gcf, 'Position', get(0,'screensize'))

savefig(gcf,fullfile(output_filepath,[bemobil_config.filename_prefix num2str(subject) '_preprocessed_interpolated_channels.fig']))
print(gcf,fullfile(output_filepath,[bemobil_config.filename_prefix num2str(subject) '_preprocessed_interpolated_channels.png']),'-dpng')
close
