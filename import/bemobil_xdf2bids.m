function bemobil_xdf2bids(config, varargin)
% Wrapper for fieldtrip function "data2bids"
% specifically for converting multimodal .xdf files to BIDS 
%
% Inputs :
%   config [struct, with required fields filename, bids_target_folder, subject, eeg.stream_keywords
%
%       config.filename               = 'P:\...SPOT_rotation\0_source-data\vp-1'\vp-1_control_body.xdf'; % required
%       config.bids_target_folder     = 'P:\...SPOT_rotation\1_BIDS-data';                               % required
%       config.subject                = 1;                                  % required
%       config.session                = 'VR';                               % optional 
%       config.run                    = 1;                                  % optional
%       config.task                   = 'rotation';                         % optional 
% 
%       config.eeg.stream_name        = 'BrainVision';                      % required
%       config.eeg.chanloc            = 'P:\...SPOT_rotation\0_raw-data\vp-1'\vp-1.elc'; % optional
%       config.eeg.new_chans          = '';                                 % optional
%         
%       config.motion.streams{1}.stream_name        = 'rigidbody1'; 
%       config.motion.streams{1}.tracking_system    = 'HTCVive'; 
%       config.motion.streams{1}.tracked_points     = 'leftFoot'; 
%       config.motion.streams{2}.stream_name        = 'rigidbody2'; 
%       config.motion.streams{2}.tracking_system    = 'HTCVive'; 
%       config.motion.streams{2}.tracked_points     = 'rightFoot'; 
%       config.motion.streams{3}.stream_name        = 'rigidbody3'; 
%       config.motion.streams{3}.tracking_system    = 'phaseSpace'; 
%       config.motion.streams{3}.tracked_points     = {'leftFoot', 'rightFoot'}; 
%
%       config.physio.streams{1}.stream_name        = {'force1'};           % optional
%
%--------------------------------------------------------------------------
% Optional Inputs :
%       Provide optional inputs as key value pairs. See below for example
%       bemobil_xdf2bids(config, 'general_metadata', generalInfo); 
% 
%       general_metadata
%       participant_metadata
%       eeg_metadata
%       motion_metadata
%       physio_metadata
% 
% Authors : 
%       Sein Jeung (seinjeung@gmail.com) & Soeren Grothkopp (email)
%--------------------------------------------------------------------------

% add load_xdf to path 
[filepath,~,~] = fileparts(which('ft_defaults'));
addpath(fullfile(filepath, 'external', 'xdf'))

%% 
%--------------------------------------------------------------------------
%                   Check import configuration 
%--------------------------------------------------------------------------

% check which modalities are included
%--------------------------------------------------------------------------
importEEG           = isfield(config, 'eeg');                               % assume EEG is always in 
importMotion        = isfield(config, 'motion');
importPhys          = isfield(config, 'phys');

if ~importEEG
    error('Importing scripts require the EEG stream to be there for event processing')
end

% check for mandatory fields 
%--------------------------------------------------------------------------
config = checkfield(config, 'filename', 'required', ''); 
config = checkfield(config, 'bids_target_folder', 'required', ''); 
config = checkfield(config, 'subject', 'required', ''); 
config.eeg = checkfield(config.eeg, 'stream_name', 'required', '');         % for now, the EEG stream has to be there for smooth processing 

% assign default values to optional fields
%--------------------------------------------------------------------------
config = checkfield(config, 'task', 'DefaultTask', 'DefaultTask');

% validate file name parts 
%--------------------------------------------------------------------------
pat = {' ' '_'};

if contains(config.task, pat)
    error('Task label MUST NOT contain space or underscore. Please change task label.')
end

if isfield(config, 'session')
    if contains(config.session, pat)
        error('Session label MUST NOT contain space or underscore. Please change task label.')
    end
end

% motion-related fields
%--------------------------------------------------------------------------
if importMotion
    
    config.motion = checkfield(config.motion, 'streams', 'required', '');
    
    for Si = 1:numel(config.motion.streams)
        config.motion.streams{Si} = checkfield(config.motion.streams{Si}, 'stream_name', 'required', '');
        config.motion.streams{Si} = checkfield(config.motion.streams{Si}, 'tracking_system', 'required', '');
        config.motion.streams{Si} = checkfield(config.motion.streams{Si}, 'tracked_points', 'required', '');
    end
    
    if isfield(config, 'bids_motionconvert_custom')
        if isempty(config.bids_motionconvert_custom)
            % funcions that resolve dataset-specific problems
            motionCustom            = 'bemobil_bids_motionconvert';
        else
            motionCustom            = config.bids_motionconvert_custom;
        end
    else
        motionCustom = 'bemobil_bids_motionconvert'; 
    end
    
end

% physio-related fields
%--------------------------------------------------------------------------
if importPhys
    
    config.phys = checkfield(config.phys, 'streams', 'required', '');
    
    for Si = 1:numel(config.phys.streams)
        config.phys.streams{Si} = checkfield(config.phys.streams{Si}, 'stream_name', 'required', '');
    end
    
    % no custom function for physio processing supported yet
    physioCustom        = 'bemobil_bids_physioconvert';

end


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

% check input from variable arguments
if ~isfield(motionInfo.motion , 'tracksys_in_session')
    warning('motionInfo.motion.tracksys_in_session not defined. Assuming all trackingsystems are present in all sessions')
    motionInfo.motion.tracksys_in_session    = true(numel(bemobil_config.session_names),numel(motionInfo.motion.trsystems));
elseif ~islogical(motionInfo.motion.tracksys_in_session)
        motionInfo.motion.tracksys_in_session = logical(motionInfo.motion.tracksys_in_session);
end

if ~isfield(motionInfo.motion, 'rb_prefix')
    error('motionInfo.motion.rb_prefix must exist')
elseif isempty(motionInfo.motion.rb_prefix)
    error('motionInfo.motion.rb_prefix must contain entries')
end

for tsi = 1:numel(motionInfo.motion.trsystems)
    if ~isfield(motionInfo.motion.TrackingSystems.(motionInfo.motion.trsystems{tsi}), 'SamplingFrequencyNominal')
        warning(['Nominal sampling frequency for ' motionInfo.motion.trsystems{tsi} ' will be defined as n/a'])
        motionInfo.motion.TrackingSystems.(motionInfo.motion.trsystems{tsi}).SamplingFrequencyNominal = 'n/a';
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
% determine number of days for shifting acq_time
shift = randi([-1000,1000]);

% create key value maps for tracksys and rb
kv_trsys_to_rb  = containers.Map(motionInfo.motion.trsystems, motionInfo.motion.rb_prefix);

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
  
        
        % loop over files in each session.
        % Here 'di' will index files as runs.
        for di = 1:numel(sortedFileNames)
            
            motionStreamNames                       = bemobil_config.rigidbody_streams; % reset after modified for multisystem use 
            
            
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
           
            % shift acq_time and look for missing acquisition time data
            if ~isfield(cfg, 'acq_times')
                warning(' acquisition times are not defined. Unable to display acq_time for eeg data.')
            elseif isempty(cfg.acq_times)
                warning(' acquisition times are not specified. Unable to display acq_time for eeg data.')
            else             
                acq_time = cfg.acq_times{pi,si};  
                acq_time([1:4]) = num2str(bemobil_config.bids_shift_acquisition_time);
                eegcfg.acq_time = datestr(datenum(acq_time) + shift,'yyyy-mm-ddTHH:MM:SS.FFF'); % microseconds are rounded
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
                        if ~isfield(motionInfo.motion , 'trsystems')
                            error('Trackingsystems must be specified. Please create motionInfo.motion.trsystems containing names of tracking systems.') 
                        elseif isempty(motionInfo.motion.trsystems)
                            error('Trackingsystems must be specified. Please enter name of tracking systems in motionInfo.motion.trsystems .') 
                        else 
                            if any(contains(motionInfo.motion.trsystems, '_'))
                                error('Name of trackingsystem MUST NOT contain underscores. Please change name of trackingsystem.')
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
                        
                        % tracksys in session
                        trsystems                       = motionInfo.motion.trsystems;
                        trsystems_in_session            = trsystems(motionInfo.motion.tracksys_in_session(si,:));  
                                              
                        % order fieldtrip data according to tracksys
                        rb_prefix_in_session = motionInfo.motion.rb_prefix(motionInfo.motion.tracksys_in_session(si,:));
                        
                        idx_ftmotion =[];
                        for iP = 1:numel(rb_prefix_in_session)
                            has_prefix = [];
                            for iM = 1:numel(ftmotion)
                                has_prefix(iM) = any(contains(ftmotion{iM}.label,rb_prefix_in_session{iP}));
                            end
                            idx_ftmotion = [idx_ftmotion find(has_prefix)];
                        end
                        
                        ftmotion = ftmotion(idx_ftmotion);
                        
                        %--------------------------------------------------
                        % determine use case
                        
                        multisys = false;
                        singlesys = false;
                        singlestream = false;
                        
                        if numel(trsystems_in_session) > 1
                            multisys = true;
                        else
                            singlesys = true;
                        end 
                        
                        if any(contains(rb_prefix_in_session,bemobil_config.rb_prefix_single_stream ))
                           singlestream = true;
                        end 
                        %--------------------------------------------------
                        % Preparing data for singlestream and multisys/singlesys case
                        if singlestream
                            
                            % checks
                            if ~isfield(bemobil_config, 'rigidbody_single_stream_names')
                                error('bemobil_config.rigidbody_single_stream_names must exist')
                            elseif isempty(bemobil_config.rigidbody_single_stream_names)
                                error('bemobil_config.rigidbody_single_stream_names must contain entries')
                            end

                            single_stream_rb_names = bemobil_config.rigidbody_single_stream_names;
                            
                            for rbi = 1:numel(bemobil_config.rb_prefix_single_stream)
                                
                                % find single stream
                                has_single_stream = [];
                                idx_single_stream = [];

                                for fti = 1:numel(ftmotion)
                                    has_single_stream(fti) = any(contains(ftmotion{fti}.label,bemobil_config.rb_prefix_single_stream{rbi}));
                                end 
                                
                                idx_single_stream = find(has_single_stream);

                                % add new motion stream names
                                motionStreamNames = motionStreamNames(bemobil_config.bids_rb_in_sessions(si,:));
                                motionStreamNames(idx_single_stream) = [];
                                
                                idx_rb_names = [];
                                idx_rb_names = find(contains(single_stream_rb_names, ...
                                    bemobil_config.rb_prefix_single_stream(rbi))); 
                                
                                % append rb names from single stream
                                motionStreamNames_in_session = [motionStreamNames single_stream_rb_names(idx_rb_names)]; 
                                
                                % check rb of single stream in session
                                if ~isfield(bemobil_config, 'bids_rb_single_stream_in_sessions')
                                    warning('bids_rb_single_stream_in_sessions not defined. Assuming all rigidbodies are present in all sessions')
                                    bemobil_config.bids_rb_single_stream_in_sessions    = true(numel(bemobil_config.session_names),numel(motionStreamNames));
                                else
                                    if ~islogical(bemobil_config.bids_rb_single_stream_in_sessions)
                                        bemobil_config.bids_rb_single_stream_in_sessions = logical(bemobil_config.bids_rb_single_stream_in_sessions);
                                    end
                                end

                                % convert single stream rbs to streams
                                ftmotion_single_stream = [];
                                ftmotion_single_stream = extractRB(ftmotion{idx_single_stream},single_stream_rb_names(idx_rb_names));

                                % remove global single stream and replace with rb as streams
                                ftmotion{idx_single_stream} = ftmotion_single_stream;

                                % group rbs in cell according to tracksys
                                indices_single_stream(rbi)= idx_single_stream;

                                for iM = 1:numel(ftmotion)
                                    if iM~=indices_single_stream
                                        ftmotion{iM} = {ftmotion{iM}};
                                    end 
                                end 
                            end
                            
                            ftmotion_grouped = ftmotion;
                                
                            % group motionStreamNames in cell according to tracksys
                            for iP = 1:numel(rb_prefix_in_session)                           
                                idx_rb = find(contains(motionStreamNames_in_session,rb_prefix_in_session{iP}));
                                motionStreamNames_grouped{iP} = motionStreamNames_in_session(idx_rb);
                            end
                            
                            motionStreamNames = motionStreamNames_grouped;
                            
                        else % multisys or singlesys without single stream
                            
                            motionStreamNames_in_session = motionStreamNames(bemobil_config.bids_rb_in_sessions(si,:));
                            % group ftmotion and motionStreamNames in cell according to tracksys
                            for iP = 1:numel(rb_prefix_in_session)
                                % ftmotion
                                has_prefix = [];
                                for iM = 1:numel(ftmotion)
                                    has_prefix(iM) = any(contains(ftmotion{iM}.label,rb_prefix_in_session{iP}));
                                end
                                idx_ftmotion = find(has_prefix);
                                ftmotion_grouped{iP} = ftmotion(idx_ftmotion);
                                
                                % motionStreamNames
                                idx_rb = find(contains(motionStreamNames_in_session,rb_prefix_in_session{iP}));
                                motionStreamNames_grouped{iP} = motionStreamNames_in_session(idx_rb);
                            end
                            
%                             ftmotion = ftmotion_grouped; 
%                             motionStreamNames = motionStreamNames_grouped;

                        end 
                        

                        %--------------------------------------------------
                            if multisys
   
                                MotionChannelCount = 0; 
                                TrackedPointsCountTotal = 0;
                          
                                for tsi = 1:numel(trsystems_in_session)
                                    
                                    tracksys          = trsystems_in_session{tsi};
                                    ftmotion          = ftmotion_grouped{tsi};
                                    motionStreamNames = motionStreamNames_grouped{tsi};

                                    % if needed, execute a custom function for any alteration to the data to address dataset specific issues
                                    % (quat2eul conversion, unwrapping of angles, resampling, wrapping back to [pi, -pi], and concatenating for instance)
                                    motion = feval(motionCustom, ftmotion, motionStreamNames, participantNr, si, di);

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

                                    motioncfg.TrackingSystemCount   = numel(trsystems_in_session);

                                    % sampling frequency
                                     if isfield(motion, 'fsample')
                                        motionInfo.motion.TrackingSystems.(tracksys).SamplingFrequencyEffective = motion.fsample;
                                     else 
                                        motionInfo.motion.TrackingSystems.(tracksys).SamplingFrequencyEffective = motion.hdr.Fs;
                                     end 

                                     if strcmpi(motionInfo.motion.TrackingSystems.(tracksys).SamplingFrequencyNominal, 'n/a')
                                        motionInfo.motion.TrackingSystems.(tracksys).SamplingFrequencyNominal = motion.hdr.nFs;
                                     end 

                                    % data type and acquisition label
                                    motioncfg.acq                                     = motionInfo.acq;

                                    % motion specific fields in json
                                    motioncfg.motion                                  = motionInfo.motion;

                                    % tracking system
                                    motioncfg.motion.trsystems                        = trsystems ; % needed for removing general trackingsys info 
                                    motioncfg.tracksys                                = tracksys; 
                                    motioncfg.motion.trsystems_in_session             = trsystems_in_session;

                                    % start time
                                    motioncfg.motion.start_time                       = motionStartTime - eegStartTime;

                                    % coordinate system
                                    motioncfg.coordsystem.MotionCoordinateSystem      = motionInfo.coordsystem.MotionCoordinateSystem;
                                    motioncfg.coordsystem.MotionRotationRule          = motionInfo.coordsystem.MotionRotationRule;
                                    motioncfg.coordsystem.MotionRotationOrder         = motionInfo.coordsystem.MotionRotationOrder;

                                    %--------------------------------------------------
                                    % rename and fill out motion-specific fields to be used in channels_tsv
                                    if singlestream
                                        rb_streams = horzcat(motionStreamNames_grouped{:}); 
                                        rb_names = bemobil_config.rigidbody_names(find(1 == bemobil_config.bids_rb_single_stream_in_sessions(si,:))); % usually rb in session would be used. Exception bc of how data is structured
                                        rb_anat = bemobil_config.rigidbody_anat(find(1 == bemobil_config.bids_rb_single_stream_in_sessions(si,:)));
                                    else
                                        rb_streams = bemobil_config.rigidbody_streams(find(1 == bemobil_config.bids_rb_in_sessions(si,:)));
                                        rb_names = bemobil_config.rigidbody_names(find(1 == bemobil_config.bids_rb_in_sessions(si,:)));
                                        rb_anat = bemobil_config.rigidbody_anat(find(1 == bemobil_config.bids_rb_in_sessions(si,:)));
                                    end


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
                                        motioncfg.channels.tracking_system{ci}     = tracksys; 


                                        % assign object names and anatomical positions

                                        for iRB = 1:numel(rb_streams)
                                            if contains(motion.hdr.label{ci},rb_streams{iRB})
                                                motioncfg.channels.tracked_point{ci}       = rb_names{iRB};
                                                if iscell(bemobil_config.rigidbody_anat)
                                                    motioncfg.channels.placement{ci}  = rb_anat{iRB};
                                                else
                                                    motioncfg.channels.placement{ci} =  rb_anat;
                                                end
                                            end

                                        end

                                        motioncfg.channels.component{ci}    = splitlabel{end}; % REQUIRED. Component of the representational system that the channel contains.            
                                    end

                                     % shift acq_time and look for missing acquisition time data
                                     if ~isfield(cfg, 'acq_times')
                                        warning(' acquisition times are not defined. Unable to display acq_time for motion data.')
                                     elseif isempty(cfg.acq_times)
                                        warning(' acquisition times are not specified. Unable to display acq_time for motion data.')
                                     else
                                        acq_time = cfg.acq_times{pi,si};
                                        acq_time([1:4]) = num2str(bemobil_config.bids_shift_acquisition_time);
                                        acq_time = datenum(acq_time) - (motioncfg.motion.start_time/(24*60*60));
                                        motioncfg.acq_time = datestr(acq_time + shift,'yyyy-mm-ddTHH:MM:SS.FFF'); % microseconds are rounded 
                                     end 

                                    % RecordingDuration
                                    fs_effective = motionInfo.motion.TrackingSystems.(tracksys).SamplingFrequencyEffective;
                                    motioncfg.motion.TrackingSystems.(tracksys).RecordingDuration = (motion.hdr.nSamples*motion.hdr.nTrials)/fs_effective;

                                    % tracked points per trackingsystem
                                    motioncfg.motion.tracksys = [];
                                    if checkequal(motionStreamNames) % checks if array contains similar entries
                                        warning(' rigidbody streams have the same name. Assuming TrackedPointsCount per trackingsystem is 1.')
                                        for ti=1:numel(motioncfg.motion.trsystems)
                                            tracksys = motioncfg.motion.trsystems{ti};
                                            motioncfg.motion.tracksys.(tracksys).TrackedPointsCount = 1; % hard coded TrackedPointsCount
                                        end 
                                    else
                                        for ti=1:numel(motioncfg.motion.trsystems)
                                            tracksys = motioncfg.motion.trsystems{ti};
                                            rb_name = kv_trsys_to_rb(tracksys); % select rigid body name corresponding to trackingsystem
                                            motioncfg.motion.tracksys.(tracksys).TrackedPointsCount = sum(contains(motionStreamNames, rb_name)); % add entries which contain rb_name for corresponding tracking system
                                        end 
                                    end 

                                    % match channel tokens with tracked points
                                    for tpi = 1:numel(bemobil_config.rigidbody_names)
                                        tokens{tpi} = ['t' num2str(tpi)];
                                    end 
                                    motioncfg.motion.tpPairs  = containers.Map(bemobil_config.rigidbody_names,tokens);

                                    % MotionChannelCount
                                    MotionChannelCount = MotionChannelCount + motion.hdr.nChans;
                                    motion.hdr.nChansTs = MotionChannelCount;

                                    % TrackedPointsCountTotal
                                    TrackedPointsCountTotal = TrackedPointsCountTotal + numel(unique(motioncfg.channels.tracked_point));
                                    motioncfg.motion.TrackedPointsCountTotal = TrackedPointsCountTotal;


                                    % write motion files in bids format
                                    data2bids(motioncfg, motion);

                                end 
                                
                            end 
                            
                            %--------------------------------------------------    
                            if singlesys
                                
                                tracksys          = trsystems_in_session{1};
                                ftmotion          = ftmotion_grouped{1};
                                motionStreamNames = motionStreamNames_grouped{1};
                                                              
                                % if needed, execute a custom function for any alteration to the data to address dataset specific issues
                                % (quat2eul conversion, unwrapping of angles, resampling, wrapping back to [pi, -pi], and concatenating for instance)
                                motion = feval(motionCustom, ftmotion, motionStreamNames, participantNr, si, di);

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

                                motioncfg.TrackingSystemCount   = numel(trsystems_in_session);

                                % sampling frequency
                                if isfield(motion, 'fsample')
                                   motionInfo.motion.TrackingSystems.(tracksys).SamplingFrequencyEffective = motion.fsample;
                                else 
                                   motionInfo.motion.TrackingSystems.(tracksys).SamplingFrequencyEffective = motion.hdr.Fs;
                                end 
                                
                                if strcmpi(motionInfo.motion.TrackingSystems.(tracksys).SamplingFrequencyNominal, 'n/a')
                                   motionInfo.motion.TrackingSystems.(tracksys).SamplingFrequencyNominal = motion.hdr.nFs;
                                end 


                                % data type and acquisition label
                                motioncfg.acq                                     = motionInfo.acq;

                                % motion specific fields in json
                                motioncfg.motion                                  = motionInfo.motion;

                                % tracking system
                                motioncfg.motion.trsystems                        = trsystems ; % needed for removing general trackingsys info 
                                motioncfg.tracksys                                = tracksys; % has to be adjusted for multiple tracking systems in one session
                                motioncfg.motion.trsystems_in_session             = trsystems_in_session;
                                
                                % start time
                                motioncfg.motion.start_time                       = motionStartTime - eegStartTime;

                                % coordinate system
                                motioncfg.coordsystem.MotionCoordinateSystem      = motionInfo.coordsystem.MotionCoordinateSystem;
                                motioncfg.coordsystem.MotionRotationRule          = motionInfo.coordsystem.MotionRotationRule;
                                motioncfg.coordsystem.MotionRotationOrder         = motionInfo.coordsystem.MotionRotationOrder;

                                %--------------------------------------------------
                                % rename and fill out motion-specific fields to be used in channels_tsv

                                if singlestream
                                    rb_streams = horzcat(motionStreamNames_grouped{:}); 
                                    rb_names = bemobil_config.rigidbody_names(find(1 == bemobil_config.bids_rb_single_stream_in_sessions(si,:))); % usually rb in session would be used. Exception bc of how data is structured
                                    rb_anat = bemobil_config.rigidbody_anat(find(1 == bemobil_config.bids_rb_single_stream_in_sessions(si,:)));
                                else
                                    rb_streams = bemobil_config.rigidbody_streams(find(1 == bemobil_config.bids_rb_in_sessions(si,:)));
                                    rb_names = bemobil_config.rigidbody_names(find(1 == bemobil_config.bids_rb_in_sessions(si,:)));
                                    rb_anat = bemobil_config.rigidbody_anat(find(1 == bemobil_config.bids_rb_in_sessions(si,:)));
                                end
                                
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
                                    motioncfg.channels.tracking_system{ci}     = motioncfg.tracksys; 


                                    % assign object names and anatomical positions
                                    for iRB = 1:numel(rb_streams)
                                        if contains(motion.hdr.label{ci}, rb_streams{iRB})
                                            motioncfg.channels.tracked_point{ci}       = rb_names{iRB};
                                            if iscell(bemobil_config.rigidbody_anat)
                                                motioncfg.channels.placement{ci}  = rb_anat{iRB};
                                            else
                                                motioncfg.channels.placement{ci} =  rb_anat;
                                            end
                                        end

                                    end

                                    motioncfg.channels.component{ci}    = splitlabel{end}; % REQUIRED. Component of the representational system that the channel contains.            
                                end



                                 % shift acq_time and look for missing acquisition time data
                                 if ~isfield(cfg, 'acq_times')
                                    warning(' acquisition times are not defined. Unable to display acq_time for motion data.')
                                 elseif isempty(cfg.acq_times)
                                    warning(' acquisition times are not specified. Unable to display acq_time for motion data.')
                                 else
                                    acq_time = cfg.acq_times{pi,si};
                                    acq_time([1:4]) = num2str(bemobil_config.bids_shift_acquisition_time);
                                    acq_time = datenum(acq_time) - (motioncfg.motion.start_time/(24*60*60));
                                    motioncfg.acq_time = datestr(acq_time + shift,'yyyy-mm-ddTHH:MM:SS.FFF'); % microseconds are rounded 
                                 end 
                                 
                                % RecordingDuration
                                fs_effective = motionInfo.motion.TrackingSystems.(tracksys).SamplingFrequencyEffective;
                                motioncfg.motion.TrackingSystems.(tracksys).RecordingDuration = (motion.hdr.nSamples*motion.hdr.nTrials)/fs_effective;
                                
                                % tracked points per trackingsystem
                                motioncfg.motion.tracksys = [];
                                if checkequal(motionStreamNames) % checks if array contains similar entries
                                    warning(' rigidbody streams have the same name. Assuming TrackedPointsCount per trackingsystem is 1.')
                                    for ti=1:numel(motioncfg.motion.trsystems)
                                        tracksys = motioncfg.motion.trsystems{ti};
                                        motioncfg.motion.tracksys.(tracksys).TrackedPointsCount = 1; % hard coded TrackedPointsCount
                                    end 
                                else
                                    for ti=1:numel(motioncfg.motion.trsystems)
                                        tracksys = motioncfg.motion.trsystems{ti};
                                        rb_name = kv_trsys_to_rb(tracksys); % select rigid body name corresponding to trackingsystem
                                        motioncfg.motion.tracksys.(tracksys).TrackedPointsCount = sum(contains(motionStreamNames, rb_name)); % add entries which contain rb_name for corresponding tracking system
                                    end 
                                end 

                                % match channel tokens with tracked points
                                for tpi = 1:numel(bemobil_config.rigidbody_names)
                                    tokens{tpi} = ['t' num2str(tpi)];
                                end 
                                motioncfg.motion.tpPairs  = containers.Map(bemobil_config.rigidbody_names,tokens);
                                
                                % MotionChannelCount
                                motion.hdr.nChansTs = motion.hdr.nChans;
                                
                                % TrackedPointCountTotal
                                motioncfg.motion.TrackedPointsCountTotal =  numel(unique(motioncfg.channels.tracked_point));
                                
                               
                                % write motion files in bids format
                                data2bids(motioncfg, motion);
                        end
                        
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


function newconfig = unit_check(oldconfig)

newconfig   = oldconfig;

% position
%-----------------------------------------------------------
    if isfield(oldconfig, 'bids_motion_position_units')
        if ~iscell(oldconfig.bids_motion_position_units)
            newconfig.bids_motion_position_units = {newconfig.bids_motion_position_units};
        end

        if numel(oldconfig.bids_motion_position_units) ~= numel(oldconfig.session_names)
            if numel(oldconfig.bids_motion_position_units) == 1
                newconfig.bids_motion_position_units = repmat(newconfig.bids_motion_position_units, 1, numel(oldconfig.session_names));
                warning('Only one pos unit specified for multiple sessions - applying same unit to all sessions')
            else
                error('Config field bids_motion_position_units must have either one entry or the number of entries (in cell array) have to match number of entries in field session_names')
            end
        end
    else
        newconfig.bids_motion_position_units       = repmat({'m'},1,numel(oldconfig.session_names));
        warning('Config field bids_motion_position_units unspecified - assuming meters')
    end
 
% orientation
%-----------------------------------------------------------
    if isfield(oldconfig, 'bids_motion_orientation_units')
        if ~iscell(oldconfig.bids_motion_orientation_units)
            newconfig.bids_motion_orientation_units = {newconfig.bids_motion_orientation_units};
        end

        if numel(oldconfig.bids_motion_orientation_units) ~= numel(oldconfig.session_names)
            if numel(oldconfig.bids_motion_orientation_units) == 1
                newconfig.bids_motion_orientation_units = repmat(newconfig.bids_motion_orientation_units, 1, numel(oldconfig.session_names));
                warning('Only one orientation unit specified for multiple sessions - applying same unit to all sessions')
            else
                error('Config field bids_motion_orientation_units must have either one entry or the number of entries (in cell array) have to match number of entries in field session_names')
            end
        end
    else
        newconfig.bids_motion_orientation_units       = repmat({'rad'},1,numel(oldconfig.session_names));
        warning('Config field bids_motion_orientation_units unspecified - assuming radians')
    end

% velocity
%-----------------------------------------------------------
    if isfield(oldconfig, 'bids_motion_velocity_units')
        if ~iscell(oldconfig.bids_motion_velocity_units)
            newconfig.bids_motion_velocity_units = {newconfig.bids_motion_velocity_units};
        end

        if numel(oldconfig.bids_motion_velocity_units) ~= numel(oldconfig.session_names)
            if numel(oldconfig.bids_motion_velocity_units) == 1
                newconfig.bids_motion_velocity_units = repmat(newconfig.bids_motion_velocity_units, 1, numel(oldconfig.session_names));
                warning('Only one orientation unit specified for multiple sessions - applying same unit to all sessions')
            else
                error('Config field bids_motion_velocity_units must have either one entry or the number of entries (in cell array) have to match number of entries in field session_names')
            end
        end
    else
        newconfig.bids_motion_velocity_units       = repmat({'m/s'},1,numel(oldconfig.session_names));
        warning('Config field bids_motion_velocity_units unspecified - assuming meters per second')
    end
    
% angularvelocity
%-----------------------------------------------------------
    if isfield(oldconfig, 'bids_motion_angularvelocity_units')
        if ~iscell(oldconfig.bids_motion_angularvelocity_units)
            newconfig.bids_motion_angularvelocity_units = {newconfig.bids_motion_angularvelocity_units};
        end

        if numel(oldconfig.bids_motion_angularvelocity_units) ~= numel(oldconfig.session_names)
            if numel(oldconfig.bids_motion_angularvelocity_units) == 1
                newconfig.bids_motion_angularvelocity_units = repmat(newconfig.bids_motion_angularvelocity_units, 1, numel(oldconfig.session_names));
                warning('Only one orientation unit specified for multiple sessions - applying same unit to all sessions')
            else
                error('Config field bids_motion_angularvelocity_units must have either one entry or the number of entries (in cell array) have to match number of entries in field session_names')
            end
        end
    else
        newconfig.bids_motion_angularvelocity_units       = repmat({'rad/s'},1,numel(oldconfig.session_names));
        warning('Config field bids_motion_angularvelocity_units unspecified - assuming radians per second')
    end
    
% acceleration
%-----------------------------------------------------------
    if isfield(oldconfig, 'bids_motion_acceleration_units')
        if ~iscell(oldconfig.bids_motion_acceleration_units)
            newconfig.bids_motion_acceleration_units = {newconfig.bids_motion_acceleration_units};
        end

        if numel(oldconfig.bids_motion_acceleration_units) ~= numel(oldconfig.session_names)
            if numel(oldconfig.bids_motion_acceleration_units) == 1
                newconfig.bids_motion_acceleration_units = repmat(newconfig.bids_motion_acceleration_units, 1, numel(oldconfig.session_names));
                warning('Only one orientation unit specified for multiple sessions - applying same unit to all sessions')
            else
                error('Config field bids_motion_acceleration_units must have either one entry or the number of entries (in cell array) have to match number of entries in field session_names')
            end
        end
    else
        newconfig.bids_motion_acceleration_units       = repmat({'m/s^2'},1,numel(oldconfig.session_names));
        warning('Config field bids_motion_acceleration_units unspecified - assuming meters per square second')
    end

% angularacceleration
%-----------------------------------------------------------
    if isfield(oldconfig, 'bids_motion_angularacceleration_units')
        if ~iscell(oldconfig.bids_motion_angularacceleration_units)
            newconfig.bids_motion_angularacceleration_units = {newconfig.bids_motion_angularacceleration_units};
        end

        if numel(oldconfig.bids_motion_angularacceleration_units) ~= numel(oldconfig.session_names)
            if numel(oldconfig.bids_motion_angularacceleration_units) == 1
                newconfig.bids_motion_angularacceleration_units = repmat(newconfig.bids_motion_angularacceleration_units, 1, numel(oldconfig.session_names));
                warning('Only one orientation unit specified for multiple sessions - applying same unit to all sessions')
            else
                error('Config field bids_motion_angularacceleration_units must have either one entry or the number of entries (in cell array) have to match number of entries in field session_names')
            end
        end
    else
        newconfig.bids_motion_angularacceleration_units       = repmat({'rad/s^2'},1,numel(oldconfig.session_names));
        warning('Config field bids_motion_angularacceleration_units unspecified - assuming radians per square second')
    end

% mangeticfield
%-----------------------------------------------------------
    if isfield(oldconfig, 'bids_motion_mangeticfield_units')
        if ~iscell(oldconfig.bids_motion_mangeticfield_units)
            newconfig.bids_motion_mangeticfield_units = {newconfig.bids_motion_mangeticfield_units};
        end

        if numel(oldconfig.bids_motion_mangeticfield_units) ~= numel(oldconfig.session_names)
            if numel(oldconfig.bids_motion_mangeticfield_units) == 1
                newconfig.bids_motion_mangeticfield_units = repmat(newconfig.bids_motion_mangeticfield_units, 1, numel(oldconfig.session_names));
                warning('Only one orientation unit specified for multiple sessions - applying same unit to all sessions')
            else
                error('Config field bids_motion_mangeticfield_units must have either one entry or the number of entries (in cell array) have to match number of entries in field session_names')
            end
        end
    else
        newconfig.bids_motion_mangeticfield_units       = repmat({'T'},1,numel(oldconfig.session_names));
        warning('Config field bids_motion_mangeticfield_units unspecified - assuming Tesla')
    end
    
% jointangle
%-----------------------------------------------------------
    if isfield(oldconfig, 'bids_motion_jointangle_units')
        if ~iscell(oldconfig.bids_motion_jointangle_units)
            newconfig.bids_motion_jointangle_units = {newconfig.bids_motion_jointangle_units};
        end

        if numel(oldconfig.bids_motion_jointangle_units) ~= numel(oldconfig.session_names)
            if numel(oldconfig.bids_motion_jointangle_units) == 1
                newconfig.bids_motion_jointangle_units = repmat(newconfig.bids_motion_jointangle_units, 1, numel(oldconfig.session_names));
                warning('Only one orientation unit specified for multiple sessions - applying same unit to all sessions')
            else
                error('Config field bids_motion_jointangle_units must have either one entry or the number of entries (in cell array) have to match number of entries in field session_names')
            end
        end
    else
        newconfig.bids_motion_jointangle_units       = repmat({'rad'},1,numel(oldconfig.session_names));
        warning('Config field bids_motion_jointangle_units unspecified - assuming radians')
    end
    
end

function y = checkequal(x)
% Input 'x' should be cell array
% Output 'y' logical value true. If any input cell array index is equal to
% another else false
% Example1:
% a{1}=[1 1 0]; a{2}=[0 0 0]; a{3}=[0 0 0];
% y = checkequal(a);
% Output is y = logical(1)
% Example2:
% a{1}=[1 1 0]; a{2}=[0 1 0]; a{3}=[0 0 0];
% y = checkequal(a);
% Output is y = logical(0)
y = false;
num = numel(x);
for i = 1:num
    for j = 1:num
        if i~=j
            if isequal(x{i},x{j})
                y = true;
                return;
            end
        end
    end
end
end