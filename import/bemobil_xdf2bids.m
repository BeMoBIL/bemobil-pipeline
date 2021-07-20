function bemobil_xdf2bids(bemobil_config, numericalIDs, varargin)

% Example script for converting .xdf recordings to BIDS
%
% Inputs :
%   bemobil_config
%   example fields
%       bemobil_config.study_folder             = 'E:\Project_BIDS\example_dataset_MWM\';
%       bemobil_config.filename_prefix          = 'sub_';
%       bemobil_config.source_data_folder       = '0_source-data\';
%       bemobil_config.bids_data_folder         = '1_BIDS-data\';
%       bemobil_config.channel_locations_filename = 'VN_E1_eloc.elc'; (or bemobil_config.elec_struct = elec)
%       bemobil_config.session_names            = {'VR' 'desktop'};
%       bemobil_config.other_data_types         = {'motion', 'physio'};
%       bemobil_config.rigidbody_streams        = {'rb_1', 'rb_2', 'rb_3', 'rb_4'}; 
%       bemobil_config.rigidbody_names          = {'playerTransform','rightHand', 'leftHand', 'Torso'};
%       bemobil_config.rigidbody_anat           = {'head','right hand','left hand','back center'};
%       bemobil_config.bids_rb_in_sessions      = [1,1,1,1;1,0,0,0];
%                                                  logicals : indicate which rbstreams are present in which sessions
%       bemobil_config.physio_streams           = {'eyetrack', 'forceplate'};
%       bemobil_config.bids_phys_in_sessions    = [1,1;1,0];
%                                                  logicals : indicate which rbstreams are present in which sessions
%       bemobil_config.bids_eeg_keyword         = {'EEG'};
%                                                  a unique keyword used to identify the eeg stream in the .xdf file
%       bemobil_config.bids_task_label          = 'VNE1';
%       bemobil_config.bids_source_zeropad      = 2;
%       bemobil_config.bids_motionconvert_custom = 'motion_customfunctionname';

%
%   numericalIDs
%       array of participant numerical IDs in the data set
%
%--------------------------------------------------------------------------
% Optional Inputs :
%
%   Note on multi-session and multi-run data sets :
%
%           Entries in bemobil_config.session_names will search through the
%           raw data directory of the participant and group together
%           .xdf files with matching keyword in the name into one session
%           If there are multiple files in one session, they will be given
%           separate 'run' numbers in the file name
%           !! IMPORTANT !!
%           The order of runs rely on incremental name sorting
%           using the "Nature Order File Sorting Tool" - Stephen Cobeldick (2021). Natural-Order Filename Sort (https://www.mathworks.com/matlabcentral/fileexchange/47434-natural-order-filename-sort), MATLAB Central File Exchange. Retrieved March 1, 2021.
%
%           for example,
%                   sub-1\sub-1_VNE1_VR_rec1.xdf
%                   sub-1\sub-1_VNE1_VR_rec2.xdf
%                   sub-1\sub-1_VNE1_desktop.xdf
%           will be organized into
%                   sub-001\ses-VR\sub-001_ses-VR_task-VNE1_run-1_eeg.bdf
%                   sub-001\ses-VR\sub-001_ses-VR_task-VNE1_run-2_eeg.bdf
%                   sub-001\ses-desktop\sub-001_ses-desktop_task-VNE1_eeg.bdf
%
% Author : Sein Jeung (seinjeung@gmail.com)
%--------------------------------------------------------------------------


% input check and default value assignment
%--------------------------------------------------------------------------
bemobil_config = checkfield(bemobil_config, 'bids_data_folder', '1_BIDS-data\', '1_BIDS-data\');
bemobil_config = checkfield(bemobil_config, 'bids_eeg_keyword','EEG', 'EEG');
bemobil_config = checkfield(bemobil_config, 'bids_task_label', 'defaulttask', 'defaulttask');
bemobil_config = checkfield(bemobil_config, 'bids_source_zeropad', 0, '0');
bemobil_config = checkfield(bemobil_config, 'other_data_types', {'motion'}, 'motion');
bemobil_config = checkfield(bemobil_config, 'bids_parsemarkers_custom', [], 'none');

% motion-related fields
%--------------------------------------------------------------------------
bemobil_config = checkfield(bemobil_config, 'rigidbody_streams', {}, 'none');
bemobil_config = checkfield(bemobil_config, 'rigidbody_names', bemobil_config.rigidbody_streams, 'values from field "rigidbody_streams"');
bemobil_config = checkfield(bemobil_config, 'rigidbody_anat', 'Undefined', 'Undefined');

if ~isfield(bemobil_config, 'bids_rb_in_sessions')
    bemobil_config.bids_rb_in_sessions    = true(numel(bemobil_config.session_names),numel(bemobil_config.rigidbody_streams));
else
    if ~islogical(bemobil_config.bids_rb_in_sessions)
        bemobil_config.bids_rb_in_sessions = logical(bemobil_config.bids_rb_in_sessions);
    end
end



if isfield(bemobil_config, 'bids_motion_position_units')
    if ~iscell(bemobil_config.bids_motion_position_units)
        bemobil_config.bids_motion_position_units = {bemobil_config.bids_motion_position_units};
    end

    if numel(bemobil_config.bids_motion_position_units) ~= numel(bemobil_config.session_names)
        if numel(bemobil_config.bids_motion_position_units) == 1
            bemobil_config.bids_motion_position_units = repmat(bemobil_config.bids_motion_position_units, 1, numel(bemobil_config.session_names));
            warning('Only one pos unit specified for multiple sessions - applying same unit to all sessions')
        else
            error('Config field bids_motion_position_units must have either one entry or the number of entries (in cell array) have to match number of entries in field session_names')
        end
    end
else
    bemobil_config.bids_motion_position_units       = repmat({'m'},1,numel(bemobil_config.session_names));
    warning('Config field bids_motion_position_units unspecified - assuming meters')
end

if isfield(bemobil_config, 'bids_motion_orientation_units')
    if ~iscell(bemobil_config.bids_motion_orientation_units)
        bemobil_config.bids_motion_orientation_units = {bemobil_config.bids_motion_orientation_units};
    end

    if numel(bemobil_config.bids_motion_orientation_units) ~= numel(bemobil_config.session_names)
        if numel(bemobil_config.bids_motion_orientation_units) == 1
            bemobil_config.bids_motion_orientation_units = repmat(bemobil_config.bids_motion_orientation_units, 1, numel(bemobil_config.session_names));
            warning('Only one orientation unit specified for multiple sessions - applying same unit to all sessions')
        else
            error('Config field bids_motion_orientation_units must have either one entry or the number of entries (in cell array) have to match number of entries in field session_names')
        end
    end
else
    bemobil_config.bids_motion_orientation_units       = repmat({'rad'},1,numel(bemobil_config.session_names));
    warning('Config field bids_motion_orientation_units unspecified - assuming radians')
end


% physio-related fields
%--------------------------------------------------------------------------
bemobil_config = checkfield(bemobil_config, 'physio_streams', {}, 'none');

if ~isfield(bemobil_config, 'bids_phys_in_sessions')
    bemobil_config.bids_phys_in_sessions    = true(numel(bemobil_config.session_names),numel(bemobil_config.physio_streams));
else
    if ~islogical(bemobil_config.bids_phys_in_sessions)
        bemobil_config.bids_phys_in_sessions = logical(bemobil_config.bids_phys_in_sessions);
    end
end

% names of the steams
motionStreamNames                       = bemobil_config.rigidbody_streams;
physioStreamNames                       = bemobil_config.physio_streams;
eegStreamName                           = {bemobil_config.bids_eeg_keyword};

if contains(bemobil_config.other_data_types, 'motion')
    if isempty(bemobil_config.bids_motionconvert_custom)
        % funcions that resolve dataset-specific problems
        motionCustom            = 'bemobil_bids_motionconvert';
    else
        motionCustom            = bemobil_config.bids_motionconvert_custom;
    end
end
% no custom function for physio processing supported yet
physioCustom        = 'bemobil_bids_physioconvert';

%--------------------------------------------------------------------------
% find optional input arguments
for iVI = 1:2:numel(varargin)
    if strcmp(varargin{iVI}, 'general_metadata')
        generalInfo         = varargin{iVI+1};
    elseif strcmp(varargin{iVI}, 'participant_metadata')
        subjectInfo         = varargin{iVI+1};
    elseif strcmp(varargin{iVI}, 'motion_metadata')
        motionInfo          = varargin{iVI+1};
    elseif strcmp(varargin{iVI}, 'eeg_metadata')
        eegInfo             = varargin{iVI+1};
    elseif strcmp(varargin{iVI}, 'physio_metadata')
        physioInfo             = varargin{iVI+1};
    else
        warning('One of the optional inputs are not valid : please see help bemobil_xdf2bids')
    end
end


% check if some participant data is already in BIDS folder
skipIndices = [];
for Pi = 1:numel(numericalIDs)
    pDir = fullfile(bemobil_config.study_folder, bemobil_config.bids_data_folder, ['sub-' num2str(numericalIDs(Pi),'%03.f')]);
    if exist(pDir, 'dir')
        disp(['Subject directory ' pDir ' exists. Skipping ... '])
        skipIndices = [skipIndices, Pi];
    end
end

if ~isempty(skipIndices)
    numericalIDs(skipIndices) = [];
end

if isempty(numericalIDs)
    disp('All participant folders were already found in BIDS directory.')
    return;
end

% check if numerical IDs match subject info, if this was specified
if exist('subjectInfo','var')

    numericalIDs            = sort(numericalIDs);
    nrColInd                = find(strcmp(subjectInfo.cols, 'nr'));
    newPInfo                = {};

    % attempt to find matching rows in subject info
    for Pi = 1:numel(numericalIDs)
        pRowInd          = find(cell2mat(subjectInfo.data(:,nrColInd)) == numericalIDs(Pi),1);
        if isempty(pRowInd)
            warning(['Participant ' num2str(numericalIDs(Pi)) ' info not given : filling with n/a'])
            emptyRow         = {numericalIDs(Pi)};
            [emptyRow{2:size(subjectInfo.data,2)}] = deal('n/a');
            newPInfo(Pi,:)   = emptyRow;
        else
            newPInfo(Pi,:)   = subjectInfo.data(pRowInd,:);
        end
    end

else
    warning('Optional input participant_metadata was not entered - participant.tsv will be omitted (NOT recommended for data sharing)')
end

%--------------------------------------------------------------------------
% add load_xdf
[filepath,~,~] = fileparts(which('ft_defaults'));
addpath(fullfile(filepath, 'external', 'xdf'))

% path to sourcedata
sourceDataPath                          = fullfile(bemobil_config.study_folder, bemobil_config.source_data_folder(1:end-1));
addpath(genpath(sourceDataPath))

% general metadata that apply to all participants
if ~exist('generalInfo', 'var')

    warning('Optional input general_metadata was not entered - using default general metadata (NOT recommended for data sharing)')

    generalInfo = [];

    % root directory (where you want your bids data to be saved)
    generalInfo.bidsroot                                = fullfile(bemobil_config.study_folder, bemobil_config.bids_data_folder);

    % required for dataset_description.json
    generalInfo.dataset_description.Name                = 'Default task';
    generalInfo.dataset_description.BIDSVersion         = 'unofficial extension';

    % optional for dataset_description.json
    generalInfo.dataset_description.License             = 'n/a';
    generalInfo.dataset_description.Authors             = 'n/a';
    generalInfo.dataset_description.Acknowledgements    = 'n/a';
    generalInfo.dataset_description.Funding             = 'n/a';
    generalInfo.dataset_description.ReferencesAndLinks  = 'n/a';
    generalInfo.dataset_description.DatasetDOI          = 'n/a';

    % general information shared across modality specific json files
    generalInfo.InstitutionName                         = 'Technische Universitaet zu Berlin';
    generalInfo.InstitutionalDepartmentName             = 'Biological Psychology and Neuroergonomics';
    generalInfo.InstitutionAddress                      = 'Strasse des 17. Juni 135, 10623, Berlin, Germany';
    generalInfo.TaskDescription                         = 'Default task generated by bemobil bidstools- no metadata present';
    generalInfo.task                                    = bemobil_config.bids_task_label;

end

cfg = generalInfo;

%--------------------------------------------------------------------------
% loop over participants
for pi = 1:numel(numericalIDs)

    participantNr   = numericalIDs(pi);
    disp(['Importing .xdf for participant ' num2str(participantNr)])

    if bemobil_config.bids_source_zeropad == 0
        participantDir  = fullfile(sourceDataPath, [bemobil_config.filename_prefix num2str(participantNr)]);
    else
        participantDir  = fullfile(sourceDataPath, [bemobil_config.filename_prefix num2str(participantNr, ['%0' num2str(bemobil_config.bids_source_zeropad) '.f'])]);
    end

    % find all .xdf files for the given session in the participant directory
    participantFiles    = dir(participantDir);
    fileNameArray       = {participantFiles.name};

    % loop over sessions
    for si = 1:numel(bemobil_config.session_names)

        sessionFiles        = participantFiles(contains(fileNameArray, '.xdf') & contains(fileNameArray, bemobil_config.session_names{si})); % compares individual files with xdf and session names

        if isempty(sessionFiles)
            warning(['No files in session ' bemobil_config.session_names{si} ' found for participant ' num2str(participantNr) ', in dir ' participantDir '. Check file path and name'])
        else
            for Fi = 1:numel(sessionFiles)
                disp(['File ' sessionFiles(Fi).name ' found.'])
            end
        end

        % sort files by natural order
        sortedFileNames     = natsortfiles({sessionFiles.name});

        % randomize dates
%         date([1:4]) = num2str(1900);
%         date = datenum(date);
%         cfg.scans.acq_time = append(datestr(date + cfg.rday));
%         
       
        
        % loop over files in each session.
        % Here 'di' will index files as runs.
        for di = 1:numel(sortedFileNames)

            % construct file and participant- and file- specific config
            % information needed to construct file paths and names
            cfg.sub                                     = num2str(participantNr,'%03.f');
            cfg.dataset                                 = fullfile(participantDir, sortedFileNames{di}); % tells loadxdf where to find dataset
            cfg.ses                                     = bemobil_config.session_names{si};
            cfg.run                                     = di;
            cfg.tracksys                                = [];


            % remove session label in uni-session case
            if numel(bemobil_config.session_names) == 1
                cfg = rmfield(cfg, 'ses');
            end

            % remove session label in uni-run case
            if numel(sortedFileNames) == 1
                cfg = rmfield(cfg, 'run');
            end

            % participant information
            if exist('subjectInfo', 'var')
                allColumns      = subjectInfo.cols;
                for iCol = 1:numel(allColumns)
                    if ~strcmp(allColumns{iCol},'nr')
                        cfg.participants.(allColumns{iCol}) = subjectInfo.data{pi, iCol};
                    end
                end
            end
            

            % load and assign streams (parts taken from xdf2fieldtrip)
            %--------------------------------------------------------------
            streams                  = load_xdf(cfg.dataset);

            % initialize an array of booleans indicating whether the streams are continuous
            iscontinuous = false(size(streams));

            names = {};

            % figure out which streams contain continuous/regular and discrete/irregular data
            for i=1:numel(streams)

                names{i}           = streams{i}.info.name;

                % if the nominal srate is non-zero, the stream is considered continuous
                if ~strcmpi(streams{i}.info.nominal_srate, '0')

                    iscontinuous(i) =  true;
                    num_samples  = numel(streams{i}.time_stamps);
                    t_begin      = streams{i}.time_stamps(1);
                    t_end        = streams{i}.time_stamps(end);
                    duration     = t_end - t_begin;

                    if ~isfield(streams{i}.info, 'effective_srate')
                        % in case effective srate field is missing, add one
                        streams{i}.info.effective_srate = (num_samples - 1) / duration;
                    elseif isempty(streams{i}.info.effective_srate)
                        % in case effective srate field value is missing, add one
                        streams{i}.info.effective_srate = (num_samples - 1) / duration;
                    end

                else
                    try
                        num_samples  = numel(streams{i}.time_stamps);
                        t_begin      = streams{i}.time_stamps(1);
                        t_end        = streams{i}.time_stamps(end);
                        duration     = t_end - t_begin;

                        % if sampling rate is higher than 20 Hz,
                        % the stream is considered continuous
                        if (num_samples - 1) / duration >= 20
                            iscontinuous(i) =  true;
                            if ~isfield(streams{i}.info, 'effective_srate')
                                % in case effective srate field is missing, add one
                                streams{i}.info.effective_srate = (num_samples - 1) / duration;
                            elseif isempty(streams{i}.info.effective_srate)
                                % in case effective srate field value is missing, add one
                                streams{i}.info.effective_srate = (num_samples - 1) / duration;
                            end
                        end
                    catch
                    end
                end

            end

            xdfeeg      = streams(contains(names,eegStreamName) & iscontinuous);
            xdfmotion   = streams(contains(names,motionStreamNames(bemobil_config.bids_rb_in_sessions(si,:))) & iscontinuous);
            xdfphysio   = streams(contains(names,physioStreamNames(bemobil_config.bids_phys_in_sessions(si,:))) & iscontinuous);
            xdfmarkers  = streams(~iscontinuous);

            %--------------------------------------------------------------
            %                  Convert EEG Data to BIDS
            %--------------------------------------------------------------
            % construct fieldtrip data
            eeg        = stream2ft(xdfeeg{1});

            % save eeg start time
            eegStartTime                = eeg.time{1}(1);

            % eeg metadata construction
            %--------------------------------------------------------------
            eegcfg                              = cfg;
            eegcfg.datatype                     = 'eeg';
            eegcfg.method                       = 'convert';

            if ~isempty(bemobil_config.channel_locations_filename)
                eegcfg.coordsystem.EEGCoordinateSystem      = 'n/a';
                eegcfg.coordsystem.EEGCoordinateUnits       = 'mm';
            end

            % try to extract information from config struct if specified
            if isfield(bemobil_config, 'ref_channel')
                eegcfg.eeg.EEGReference                 = bemobil_config.ref_channel;
            end

            if isfield(bemobil_config, 'linefreqs')
                if numel(bemobil_config.linefreqs) == 1
                    eegcfg.eeg.PowerLineFrequency           = bemobil_config.linefreqs;
                elseif numel(bemobil_config.linefreqs) > 1
                    eegcfg.eeg.PowerLineFrequency           = bemobil_config.linefreqs(1);
                    warning('Only the first value specified in bemobil_config.linefreqs entered in eeg.json')
                end
            end

            % overwrite some fields if specified
            if exist('eegInfo','var')
                if isfield(eegInfo, 'eeg')
                    eegcfg.eeg          = eegInfo.eeg;
                end
                if isfield(eegInfo, 'coordsystem')
                    eegcfg.coordsystem  = eegInfo.coordsystem;
                end
            end

            % read in the event stream (synched to the EEG stream)
            if ~isempty(xdfmarkers)

                if any(cellfun(@(x) ~isempty(x.time_series), xdfmarkers))

                    events                  = stream2events(xdfmarkers, xdfeeg{1}.time_stamps);
                    eventsFound             = 1;

                    % event parser script
                    if isempty(bemobil_config.bids_parsemarkers_custom)
                        [events, eventsJSON] = bemobil_bids_parsemarkers(events);
                    else
                        [events, eventsJSON] = feval(bemobil_config.bids_parsemarkers_custom, events);
                    end

                    eegcfg.events = events;

                end
            end

            if isfield(bemobil_config, 'elec_struct')
                eegcfg.elec                         = bemobil_config.elec_struct;
            elseif isfield(bemobil_config, 'channel_locations_filename')
                eegcfg.elec = fullfile(participantDir, [bemobil_config.filename_prefix, num2str(participantNr) '_' bemobil_config.channel_locations_filename]);
            end

            % write eeg files in bids format
            data2bids(eegcfg, eeg);

            for Ti = 1:numel(bemobil_config.other_data_types)

                switch bemobil_config.other_data_types{Ti}

                    case 'motion'
                        %--------------------------------------------------
                        %            Convert Motion Data to BIDS
                        %--------------------------------------------------
                        
                        % check for missing or wrong trackingsystem information
                        if ~isfield(motionInfo.motion , 'tracksys')
                            error ('Trackingsystems must be specified. Please create motionInfo.motion.tracksys containing names of tracking systems.') 
                        elseif isempty (motionInfo.motion.tracksys)
                            error ('Trackingsystems must be specified. Please enter name of tracking systems in motionInfo.motion.tracksys .') 
                        else 
                            if any(contains(motionInfo.motion.tracksys, '_'))
                                error ('Name of trackingsystem is must not contain underscores. Please change name of trackingsystem.')
                            end
                        end
                        
                        % check if any motion data was found at all
                        if isempty(xdfmotion)
                            continue;
                        end

                        ftmotion = {};

                        % construct fieldtrip data
                        for iM = 1:numel(xdfmotion)
                            ftmotion{iM} = stream2ft(xdfmotion{iM});
                        end

                        % if needed, execute a custom function for any alteration to the data to address dataset specific issues
                        % (quat2eul conversion, unwrapping of angles, resampling, wrapping back to [pi, -pi], and concatenating for instance)
                        motion = feval(motionCustom, ftmotion{iM}, motionStreamNames(bemobil_config.bids_rb_in_sessions(si,:)), participantNr, si, di);

                        % save motion start time
                        motionStartTime              = motion.time{1}(1);

                        % construct motion metadata
                        % copy general fields
                        motioncfg       = cfg;
                        motioncfg.datatype                                = 'motion';

                        %--------------------------------------------------
                        if ~exist('motionInfo', 'var')

                            % data type and acquisition label
                            motionInfo.acq                                     = 'Motion';

                            % motion specific fields in json
                            motionInfo.motion.Manufacturer                     = 'Undefined';
                            motionInfo.motion.ManufacturersModelName           = 'Undefined';
                            motionInfo.motion.RecordingType                    = 'continuous';

                            % coordinate system
                            motionInfo.coordsystem.MotionCoordinateSystem      = 'Undefined';
                            motionInfo.coordsystem.MotionRotationRule          = 'Undefined';
                            motionInfo.coordsystem.MotionRotationOrder         = 'Undefined';

                        end

                        % sampling frequency
                        motionInfo.motion.TrackingSystems.(motionInfo.motion.tracksys{si}).SamplingFrequencyEffective = motion.hdr.Fs;
                        
                        if strcmpi(motionInfo.motion.TrackingSystems.(motionInfo.motion.tracksys{si}).SamplingFrequencyNominal, 'n/a')
                           motionInfo.motion.TrackingSystems.(motionInfo.motion.tracksys{si}).SamplingFrequencyNominal = motion.hdr.nFs;
                        end 
                        
                        
                          
                        
                        % data type and acquisition label
                        motioncfg.acq                                     = motionInfo.acq;

                        % motion specific fields in json
                        motioncfg.motion                                  = motionInfo.motion;
                        
                        % tracking system
                        motioncfg.tracksys                                = motionInfo.motion.tracksys {si}; % has to be adjusted for multiple tracking systems in one session
                        motioncfg.motion.tracksys_all                     = motionInfo.motion.tracksys ; % needed for removing general trackingsys info 
            
                        % number of all tracking systems used in session
                        motioncfg.motion.TrackingSystemCount              = numel(xdfmotion);

                        % start time
                        motioncfg.motion.start_time                       = motionStartTime - eegStartTime;

                        % coordinate system
                        motioncfg.coordsystem.MotionCoordinateSystem      = motionInfo.coordsystem;
                        

                        %--------------------------------------------------
                        % rename and fill out motion-specific fields to be used in channels_tsv
                        motioncfg.channels.name                 = cell(motion.hdr.nChans,1);
                        motioncfg.channels.tracked_point        = cell(motion.hdr.nChans,1);
                        motioncfg.channels.component            = cell(motion.hdr.nChans,1);
                        motioncfg.channels.placement            = cell(motion.hdr.nChans,1);
                        motioncfg.channels.datafile             = cell(motion.hdr.nChans,1);

                        for ci  = 1:motion.hdr.nChans

                            if  contains(motion.hdr.chantype{ci},'position')
                                motion.hdr.chantype{ci} = 'POS';
                                motion.hdr.chanunit{ci} = bemobil_config.bids_motion_position_units{si};
                            end

                            if  contains(motion.hdr.chantype{ci},'orientation')
                                motion.hdr.chantype{ci} = 'ORNT';
                                motion.hdr.chanunit{ci} = bemobil_config.bids_motion_orientation_units{si};
                            end
                            
                            if  contains(motion.hdr.chantype{ci},'velocity')
                                motion.hdr.chantype{ci} = 'VEL';
                                motion.hdr.chanunit{ci} = bemobil_config.bids_motion_orientation_units{si};
                            end
                            
                            if  contains(motion.hdr.chantype{ci},'angularvelocity')
                                motion.hdr.chantype{ci} = 'ANGVEL';
                                motion.hdr.chanunit{ci} = bemobil_config.bids_motion_orientation_units{si};
                            end
                            
                            if  contains(motion.hdr.chantype{ci},'acceleration')
                                motion.hdr.chantype{ci} = 'ACC';
                                motion.hdr.chanunit{ci} = bemobil_config.bids_motion_orientation_units{si};
                            end
                            
                            if  contains(motion.hdr.chantype{ci},'angularacceleration')
                                motion.hdr.chantype{ci} = 'ANGACC';
                                motion.hdr.chanunit{ci} = bemobil_config.bids_motion_orientation_units{si};
                            end
                            
                            if  contains(motion.hdr.chantype{ci},'magneticfield')
                                motion.hdr.chantype{ci} = 'MAGN';
                                motion.hdr.chanunit{ci} = bemobil_config.bids_motion_orientation_units{si};
                            end
                            
                            if  contains(motion.hdr.chantype{ci},'jointangle')
                                motion.hdr.chantype{ci} = 'JNTANG';
                                motion.hdr.chanunit{ci} = bemobil_config.bids_motion_orientation_units{si};
                            end

                            splitlabel                          = regexp(motion.hdr.label{ci}, '_', 'split');
                            motioncfg.channels.name{ci}         = motion.hdr.label{ci};

                            % assign object names and anatomical positions
                            for iRB = 1:numel(bemobil_config.rigidbody_streams)
                                if contains(motion.hdr.label{ci}, bemobil_config.rigidbody_streams{iRB})
                                    motioncfg.channels.tracking_system{ci}      = motionInfo.motion.tracksys {iRB};
                                    motioncfg.channels.tracked_point{ci}       = bemobil_config.rigidbody_names{iRB};
                                    if iscell(bemobil_config.rigidbody_anat)
                                        motioncfg.channels.placement{ci}  = bemobil_config.rigidbody_anat{iRB};
                                    else
                                        motioncfg.channels.placement{ci} =  bemobil_config.rigidbody_anat;
                                    end
                                end
                            end

                            motioncfg.channels.component{ci}    = splitlabel{end}; % REQUIRED. Component of the representational system that the channel contains.

                        end

                        % write motion files in bids format
                        data2bids(motioncfg, motion);

                    case 'physio'
                        %--------------------------------------------------
                        %         Convert Generic Physio Data to BIDS
                        %--------------------------------------------------
                        % check if any motion data was found at all
                        if isempty(xdfphysio)
                            continue;
                        end

                        ftphysio = {};

                        % construct fieldtrip data
                        for iP = 1:numel(xdfphysio)
                            ftphysio{iP} = stream2ft(xdfphysio{iP});
                        end

                        % resample data to match the stream of highest srate (no custom processing supported for physio data yet)
                        physio = feval(physioCustom, ftphysio, physioStreamNames(bemobil_config.bids_phys_in_sessions(si,:)), participantNr, si, di);

                        % save motion start time
                        physioStartTime              = physio.time{1}(1);

                        % construct motion metadata
                        % copy general fields
                        physiocfg               = cfg;
                        physiocfg.datatype      = 'physio';

                        %--------------------------------------------------------------
                        if ~exist('physioInfo', 'var')

                            % motion specific fields in json
                            physioInfo.physio.Manufacturer                     = 'Undefined';
                            physioInfo.physio.ManufacturersModelName           = 'Undefined';
                            physioInfo.physio.RecordingType                    = 'continuous';

                        end

                        % physio specific fields in json
                        physiocfg.physio                                  = physioInfo.physio;

                        % start time
                        physiocfg.physio.StartTime                        = physioStartTime - eegStartTime;

                        % write motion files in bids format
                        data2bids(physiocfg, physio);

                    otherwise
                        warning(['Unknown data type' bemobil_config.other_data_types{Ti}])
                end
            end
        end
    end
end







% add general json files
%--------------------------------------------------------------------------
ft_hastoolbox('jsonlab', 1);

if exist('subjectInfo', 'var')
    % participant.json
    pJSONName       = fullfile(cfg.bidsroot, 'participants.json');
    pfid            = fopen(pJSONName, 'wt');
    pString         = savejson('', subjectInfo.fields, 'NaN', '"n/a"', 'ParseLogical', true);
    fwrite(pfid, pString); fclose(pfid);
end

if eventsFound
    % events.json
    eJSONName       = fullfile(cfg.bidsroot, ['task-' cfg.task '_events.json']);
    efid            = fopen(eJSONName, 'wt');
    eString         = savejson('', eventsJSON, 'NaN', '"n/a"', 'ParseLogical', true);
    fwrite(efid, eString); fclose(efid);
end

end


%--------------------------------------------------------------------------
function [newconfig] =  checkfield(oldconfig, fieldName, defaultValue, defaultValueText)

newconfig   = oldconfig;

if ~isfield(oldconfig, fieldName)
    newconfig.(fieldName) = defaultValue;
    warning(['Config field ' fieldName ' not specified- using default value: ' defaultValueText])
end

end

%--------------------------------------------------------------------------
function [ftdata] = stream2ft(xdfstream)

% construct header
hdr.Fs                  = 'n/a';
hdr.Fs                  = xdfstream.info.effective_srate;
hdr.nFs                 = 'n/a';
hdr.nFs                 = str2num(xdfstream.info.nominal_srate);
hdr.nSamplesPre         = 0;
hdr.nSamples            = length(xdfstream.time_stamps);
hdr.nTrials             = 1;
hdr.FirstTimeStamp      = xdfstream.time_stamps(1);
hdr.TimeStampPerSample  = (xdfstream.time_stamps(end)-xdfstream.time_stamps(1)) / (length(xdfstream.time_stamps) - 1);
if isfield(xdfstream.info.desc, 'channels')
    hdr.nChans    = numel(xdfstream.info.desc.channels.channel);
else
    hdr.nChans    = str2double(xdfstream.info.channel_count);
end

hdr.label       = cell(hdr.nChans, 1);
hdr.chantype    = cell(hdr.nChans, 1);
hdr.chanunit    = cell(hdr.nChans, 1);

prefix = xdfstream.info.name;
for j=1:hdr.nChans
    if isfield(xdfstream.info.desc, 'channels')
        hdr.label{j} = [prefix '_' xdfstream.info.desc.channels.channel{j}.label];
        hdr.chantype{j} = xdfstream.info.desc.channels.channel{j}.type;
        try
            hdr.chanunit{j} = xdfstream.info.desc.channels.channel{j}.unit;
        catch
            disp([hdr.label{j} ' missing unit'])
        end
    else
        % the stream does not contain continuously sampled data
        hdr.label{j} = num2str(j);
        hdr.chantype{j} = 'unknown';
        hdr.chanunit{j} = 'unknown';
    end
end

% keep the original header details
hdr.orig = xdfstream.info;

ftdata.trial    = {xdfstream.time_series};
ftdata.time     = {xdfstream.time_stamps};
ftdata.hdr = hdr;
ftdata.label = hdr.label;

end

function outEvents = stream2events(inStreams, dataTimes)

outEvents = [];

for Si = 1:numel(inStreams)
    if iscell(inStreams{Si}.time_series)
        eventsInStream              = cell2struct(inStreams{Si}.time_series, 'value')';
        [eventsInStream.type]       = deal(inStreams{Si}.info.type);
        times                       = num2cell(inStreams{Si}.time_stamps);
        [eventsInStream.timestamp]  = times{:};
        samples                     = cellfun(@(x) find(dataTimes >= x, 1,'first'), times, 'UniformOutput', false);
        [eventsInStream.sample]     = samples{:};
        [eventsInStream.offset]     = deal([]);
        [eventsInStream.duration]   = deal([]);
        outEvents = [outEvents eventsInStream];
    end
end

% sort events by sample
[~,I] = sort([outEvents.timestamp]);
outEvents   = outEvents(I);

% re-order fields to match ft events output
outEvents   = orderfields(outEvents, [2,1,3,4,5,6]);

end
