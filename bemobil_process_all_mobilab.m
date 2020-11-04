% bemobil_process_all_mobilab - wrapper function that incorporates all necessary processing steps from raw .xdf to .set
% files in EEGLAB. Data is being loaded into mobilab, rigidbody mocap streams are processed (filtered, transformed into
% euler angles, derived) and data and marker streams are then exported to EEGLAB containing EEG and all other kinds of
% data. The dataset has the suffix '_MoBI'. This dataset is then split into individual channel types (e.g. 'EEG',
% 'MOCAP', 'EYE', 'OTHER'), and subsequently all EEG files (from several raw .xdf files) will be merged into one large
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

input_filepath = [bemobil_config.study_folder bemobil_config.raw_data_folder bemobil_config.filename_prefix num2str(subject)];
output_filepath_mobi = [bemobil_config.study_folder bemobil_config.mobilab_data_folder bemobil_config.filename_prefix num2str(subject)];
output_filepath = [bemobil_config.study_folder bemobil_config.raw_EEGLAB_data_folder bemobil_config.filename_prefix num2str(subject)];

% get rid of memory mapped object storage
try pop_editoptions( 'option_storedisk', 1, 'option_savetwofiles', 1, 'option_saveversion6', 0, 'option_single', 1, 'option_memmapdata', 0, 'option_eegobject', 0, 'option_computeica', 1, 'option_scaleicarms', 1, 'option_rememberfolder', 1, 'option_donotusetoolboxes', 0, 'option_checkversion', 1, 'option_chat', 1); catch warning('Could NOT edit EEGLAB memory options!!'); end

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

	% check if at least the mobilab part was done already
	try

		EEG = pop_loadset('filename', [bemobil_config.filename_prefix num2str(subject) '_'...
			bemobil_config.merged_filename], 'filepath', output_filepath);
		[ALLEEG, EEG_merged, CURRENTSET] = pop_newset(ALLEEG, EEG, 0,'study',0);

		warning('Old merged file already existed, using that file!')

	end
end
%% run MoBILAB

if ~exist('EEG_merged','var')
	for i_filename = 1:length(bemobil_config.filenames)
		
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
					
					% process mocap channels
					
					% quaternion values can flip their signs, which does indicate the same orientation. we need them to be
					% a smooth curve for filtering, though
					unflipped = mobilab.allStreams().item{mobilab_rb_index}.unflipSigns();
					
					% lowpass filter with the respective frequency (e.g. 6 Hz), zero-lag FIR filter
					filtered = unflipped.lowpass(bemobil_config.mocap_lowpass);
					
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
		
		%% export MoBI dataset, split into unique channel types
		
		disp('Exporting MoBI dataset. This may take a while!')
		exported_EEG = mobilab.allStreams().export2eeglab(all_data_stream_indices,all_event_stream_indices);
		
		% clear RAM
		mobilab.allStreams = [];
		
		mkdir(output_filepath)
		pop_saveset( exported_EEG, 'filename',[full_filename '_MoBI'],'filepath', output_filepath);
		disp('...done');
		clear exported_EEG
		
	end
	
	% get rid of memory mapped object storage
	try pop_editoptions( 'option_storedisk', 1, 'option_savetwofiles', 1, 'option_saveversion6', 0, 'option_single', 1, 'option_memmapdata', 0, 'option_eegobject', 0, 'option_computeica', 1, 'option_scaleicarms', 1, 'option_rememberfolder', 1, 'option_donotusetoolboxes', 0, 'option_checkversion', 1, 'option_chat', 1); catch warning('Could NOT edit EEGLAB memory options!!'); end
		
	for i_filename = 1:length(bemobil_config.filenames)
		
		full_filename = [bemobil_config.filename_prefix num2str(subject) '_' bemobil_config.filenames{i_filename}];
		
		MoBI_EEG = pop_loadset('filename',[full_filename '_MoBI.set'],'filepath', output_filepath);
		
		% split MoBI dataset into unique channel types
		[~, ~, ~, EEG_split_sets] = bemobil_split_MoBI_set(ALLEEG, MoBI_EEG, CURRENTSET);
		
		for i_splitset = 1:length(EEG_split_sets)
			pop_saveset( EEG_split_sets(i_splitset), 'filename',[full_filename '_'...
				EEG_split_sets(i_splitset).chanlocs(1).type],'filepath', output_filepath);
			disp('...done');
		end
		
		% clear RAM
		STUDY = []; CURRENTSTUDY = 0; ALLEEG = [];  CURRENTSET=[]; EEG_split_sets=[];
		
	end
	
	% This merges all EEG data files into one big file
	input_filepath = [bemobil_config.study_folder bemobil_config.raw_EEGLAB_data_folder bemobil_config.filename_prefix num2str(subject)];
	output_filepath = [bemobil_config.study_folder bemobil_config.raw_EEGLAB_data_folder bemobil_config.filename_prefix num2str(subject)];
	
	% make sure EEGLAB has no files other than the ones to be merged
	STUDY = []; CURRENTSTUDY = 0; ALLEEG = []; EEG=[]; CURRENTSET=[];
	
	% get rid of memory mapped object storage
	try pop_editoptions( 'option_storedisk', 1, 'option_savetwofiles', 1, 'option_saveversion6', 0, 'option_single', 1, 'option_memmapdata', 0, 'option_eegobject', 0, 'option_computeica', 1, 'option_scaleicarms', 1, 'option_rememberfolder', 1, 'option_donotusetoolboxes', 0, 'option_checkversion', 1, 'option_chat', 1); catch warning('Could NOT edit EEGLAB memory options!!'); end
	
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
end

%% preprocess

output_filepath = [bemobil_config.study_folder bemobil_config.raw_EEGLAB_data_folder bemobil_config.filename_prefix num2str(subject)];
if ~isempty(bemobil_config.channel_locations_filename)
	channel_locations_filepath = [bemobil_config.study_folder bemobil_config.raw_data_folder...
		bemobil_config.filename_prefix num2str(subject) '\' bemobil_config.filename_prefix num2str(subject) '_'...
		bemobil_config.channel_locations_filename];
else
	channel_locations_filepath = [];
end

% preprocessing: enter chanlocs, remove unused channels, declare EOG, resample
[ALLEEG, EEG_preprocessed, CURRENTSET] = bemobil_preprocess(ALLEEG, EEG_merged, CURRENTSET, channel_locations_filepath,...
	bemobil_config.channels_to_remove, bemobil_config.eog_channels, bemobil_config.resample_freq,...
	[bemobil_config.filename_prefix num2str(subject) '_' bemobil_config.preprocessed_filename], output_filepath,...
	bemobil_config.rename_channels, bemobil_config.ref_channel);

disp('Preprocessing done!')

% save RAM
clear EEG_merged

%% find bad channels

output_filepath = [bemobil_config.study_folder bemobil_config.raw_EEGLAB_data_folder bemobil_config.filename_prefix num2str(subject)];

% filter (necessary for findNoisyChannels)
disp('Filtering for automatic bad channel detection...')
lowcutoff = 1;
highcutoff = [];
[ ALLEEG, filteredEEG, CURRENTSET ] = bemobil_filter(ALLEEG, EEG_preprocessed, CURRENTSET, lowcutoff, highcutoff);

% Add path to prepPipeline subdirectories if not in the list
tmp = which('getPipelineDefaults');
if isempty(tmp)
	myPath = fileparts(which('prepPipeline'));
	addpath(genpath(myPath));
end

% use prep pipeline subfunction to determine noisy channels
disp('Detecting bad channels...')

noisyOut = findNoisyChannels(filteredEEG);
chans_to_interp = noisyOut.noisyChannels.all;
disp({EEG_preprocessed.chanlocs(chans_to_interp).labels})

% take EOG and REF out of the channels to interpolate: EOG is very likely
% to be different from the others but rightfully so, and the REF channel is
% currently empty, since it was added new with zeros to be filled after
% rereferencing
disp('Ignoring EOG channels for interpolation:')

disp({EEG_preprocessed.chanlocs(chans_to_interp(strcmp({EEG_preprocessed.chanlocs(chans_to_interp).type},'EOG'))).labels})
chans_to_interp(strcmp({EEG_preprocessed.chanlocs(chans_to_interp).type},'EOG'))=[];

disp('Ignoring REF channel for interpolation:')

disp({EEG_preprocessed.chanlocs(chans_to_interp(strcmp({EEG_preprocessed.chanlocs(chans_to_interp).type},'REF'))).labels})
chans_to_interp(strcmp({EEG_preprocessed.chanlocs(chans_to_interp).type},'REF'))=[];

%% do the interpolation and average referencing (reference is not considering EOGs)

disp('Interpolating bad channels...')
[ALLEEG, EEG_interp_avRef, CURRENTSET] = bemobil_interp_avref( EEG_preprocessed , ALLEEG, CURRENTSET, chans_to_interp,...
	[bemobil_config.filename_prefix num2str(subject) '_' bemobil_config.interpolated_avRef_filename], output_filepath);

