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
%       config.acquisition_time       = [2021,9,30,18,14,0.00];             % optional ([YYYY,MM,DD,HH,MM,SS]) 
% 
%       config.eeg.stream_name        = 'BrainVision';                      % required
%       config.eeg.chanloc            = 'P:\...SPOT_rotation\0_raw-data\vp-1'\vp-1.elc'; % optional
%       config.eeg.elec_struct        = elecStruct                          % optional, alternative to config.eeg.chanloc. Output struct of ft_read_sens 
%       config.eeg.chanloc_newname    = {'chan1', 'chan2'}                  % optional, cell array of size nchan X 1,  containing new chanloc labels in case you want to rename them 
%         
%       config.motion.streams{1}.stream_name        = 'stream1';            % keyword in stream name  
%                                                                                   % searched for in field "xdfdata{streamIndex}.info.name"
%       config.motion.streams{1}.tracking_system    = 'HTCVive';            % user-defined name of the tracking system
%                                                                                   % in case motion metadata are provided, match with fieldname in "motionInfo.motion.TrackingSystems.(fieldname)"
%                                                                                   % e.g., motionInfo.motion.TrackingSystems.HTCVive.Manufacturer = 'HTC'; 
%       config.motion.streams{1}.tracked_points     = 'rigid1';             % keyword in channel names, indicating which object (tracked point) is included in the stream 
%                                                                                   % searched for in field "xdfdata{streamIndex}.info.desc.channels.channel{channelIndex}.label"
%                                                                                   % required to be unique in a single tracking system
%       config.motion.streams{1}.tracked_points_anat = 'leftFoot';          % user-defined anatomical name of the point being tracked 
%       config.motion.streams{2}.stream_name        = 'stream2'; 
%       config.motion.streams{2}.tracking_system    = 'HTCVive'; 
%       config.motion.streams{2}.tracked_points     = 'rigid2';
%       config.motion.streams{2}.tracked_points_anat = 'rightFoot';
%       config.motion.streams{3}.stream_name        = 'stream3'; 
%       config.motion.streams{3}.tracking_system    = 'phaseSpace'; 
%       config.motion.streams{3}.tracked_points     = {'leftFoot', 'rightFoot'}; % keywords in channel names, when multiple tracked points are included in a single stream 
%       config.motion.POS.unit                      = 'vm';                 % optional, in case you want to use custom unit
%
%       config.phys.streams{1}.stream_name        = {'force1'};             % optional
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
%       Sein Jeung (seinjeung@gmail.com) & Soeren Grothkopp (s.grothkopp@secure.mailbox.org)
%--------------------------------------------------------------------------

% add load_xdf to path 
ft_defaults
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

% acquisition time 
config = checkfield(config, 'acquisition_time', [1800,12,31,5,5,5.000], 'default time');         % for now, the EEG stream has to be there for smooth processing 

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



% EEG-related fields
%--------------------------------------------------------------------------
if importEEG
    config.eeg = checkfield(config.eeg, 'stream_name', 'required', ''); 
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
    
    if isfield(config.motion, 'custom_function')
        if isempty(config.motion.custom_function)
            % functions that resolve dataset-specific problems
            motionCustom            = 'bemobil_bids_motionconvert';
        else
            motionCustom            = config.motion.custom_function;
        end

    else
        motionCustom = 'bemobil_bids_motionconvert'; 
    end
    
    % default channel types and units
    motion_type.POS.unit        = 'm'; 
    motion_type.ORNT.unit       = 'rad';
    motion_type.VEL.unit        = 'm/s';
    motion_type.ANGVEL.unit     = 'r/s';
    motion_type.ACC.unit        = 'm/s^2';
    motion_type.ANGACC.unit     = 'r/s^2';
    motion_type.MAGN.unit       = 'fm';
    motion_type.JNTANG.unit     = 'r';

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


%% 
%--------------------------------------------------------------------------
%                           Check metadata 
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

% check general metadata 
%--------------------------------------------------------------------------
if ~exist('generalInfo', 'var')
    
    warning('Optional input general_metadata was not entered - using default general metadata (NOT recommended for data sharing)')
    
    generalInfo = [];
    
    % root directory (where you want your bids data to be saved)
    generalInfo.bidsroot                                = config.bids_target_folder;
    
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
    generalInfo.task                                    = 'DefaultTask';
    
end

cfg = generalInfo;


% check motion config and metadata
%--------------------------------------------------------------------------
if importMotion
    
    % check how many different tracking systems are specified
    for Si = 1:numel(config.motion.streams)
        streamNames{Si}     = config.motion.streams{Si}.stream_name;
        trackSysNames{Si}   = config.motion.streams{Si}.tracking_system;
        trackedPointNames{Si} = config.motion.streams{Si}.tracked_points; 
        
        if isfield(config.motion.streams{Si}, 'tracked_points_anat')
            anatomicalNames{Si} = config.motion.streams{Si}.tracked_points_anat;
        else
            anatomicalNames{Si} = config.motion.streams{Si}.tracked_points;
        end
    end
    
    % get the (unique) tracking system names included in data
    trackSysInData   = unique(trackSysNames);
  
    % get all stream names corresponding to each tracking system
    for Si = 1:numel(trackSysInData)
        trackSysInds = find(strcmp(trackSysNames, trackSysInData{Si})); 
        streamsInData{Si}           = streamNames(trackSysInds);
        trackedPointsInData{Si}     = trackedPointNames(trackSysInds); 
        anatInData{Si}              = anatomicalNames(trackSysInds); 
    end
    
    % construct default info for tracking systems
    tracking_systems                                    = trackSysInData;
    
    for Ti = 1:numel(tracking_systems)
        defaultTrackingSystems.(tracking_systems{Ti}).Manufacturer                     = 'DefaultManufacturer';
        defaultTrackingSystems.(tracking_systems{Ti}).ManufacturersModelName           = 'DefaultModel';
        defaultTrackingSystems.(tracking_systems{Ti}).SamplingFrequencyNominal         = 'n/a'; %  If no nominal Fs exists, n/a entry returns 'n/a'. If it exists, n/a entry returns nominal Fs from motion stream.
    end
    
    if ~exist('motionInfo', 'var')
        
        warning('Optional input motion_metadata was not entered - using default metadata (NOT recommended for data sharing)')
        
        % motion specific fields in json
        motionInfo.motion = [];
        motionInfo.motion.RecordingType                     = 'continuous';
        
        % default tracking system information
        motionInfo.motion.TrackingSystems = defaultTrackingSystems;
        
        % coordinate system
        motionInfo.coordsystem.MotionCoordinateSystem      = 'RUF';
        motionInfo.coordsystem.MotionRotationRule          = 'left-hand';
        motionInfo.coordsystem.MotionRotationOrder         = 'ZXY';
        
    else
        if isfield(motionInfo, 'motion')
            if isfield(motionInfo.motion, 'TrackingSystems')
                
                % take all tracking systems defined in the metadata input
                trackSysInMeta = fieldnames(motionInfo.motion.TrackingSystems);
                
                % identify tracking systems in the data but not in metadata
                trackSysNoMeta  = setdiff(trackSysInData, trackSysInMeta);
                
                % construct metadata for ones that are missing them
                for Ti = 1:numel(trackSysNoMeta)
                    motionInfo.motion.TrackingSystems.(trackSysNoMeta{Ti}).Manufacturer                     = 'DefaultManufacturer';
                    motionInfo.motion.TrackingSystems.(trackSysNoMeta{Ti}).ManufacturerModelName            = 'DefaultModel';
                    motionInfo.motion.TrackingSystems.(trackSysNoMeta{Ti}).SamplingFrequencyNominal         = 'n/a';
                end
                
                % identify tracking systems in metadata but not in the data
                trackSysNoData = setdiff(trackSysInMeta, trackSysInData);
                
                % remove unused tracking systems from metadata struct
                motionInfo.motion.TrackingSystems = rmfield(motionInfo.motion.TrackingSystems, trackSysNoData);
            else
                warning('No information on tracking system given - filling with default info')
                
                % default tracking system information
                motionInfo.motion.TrackingSystems = defaultTrackingSystems;
            end
        else
            warning('No information for motion json given - filling with default info')
            
            % motion specific fields in json
            motionInfo.motion.RecordingType                     = 'continuous';
            
            % default tracking system information
            motionInfo.motion.TrackingSystems = defaultTrackingSystems;
        end        
    end
    
    % create key value maps for tracking systems and stream names
    kv_trsys_to_st          = containers.Map(trackSysInData, streamsInData);
    kv_trsys_to_trp         = containers.Map(trackSysInData, trackedPointsInData); 
    kv_trsys_to_anat        = containers.Map(trackSysInData, anatInData); 
end


%%
% check if numerical IDs match subject info, if this was specified
%--------------------------------------------------------------------------
if exist('subjectInfo','var')
    
    nrColInd                = find(strcmp(subjectInfo.cols, 'nr'));
    
    % attempt to find matching rows in subject info
    pRowInd          = find(cell2mat(subjectInfo.data(:,nrColInd)) == config.subject,1);
    if isempty(pRowInd)
        warning(['Participant ' num2str(numericalIDs(Pi)) ' info not given : filling with n/a'])
        emptyRow         = {config.subject};
        [emptyRow{2:size(subjectInfo.data,2)}] = deal('n/a');
        newPInfo   = emptyRow;
    else
        newPInfo   = subjectInfo.data(pRowInd,:);
    end
    
else
    warning('Optional input participant_metadata was not entered - participant.tsv will be omitted (NOT recommended for data sharing)')
end

% construct file and participant- and file- specific config
% information needed to construct file paths and names
%--------------------------------------------------------------------------
cfg.sub                                     = num2str(config.subject);
cfg.dataset                                 = config.filename; 
cfg.bidsroot                                = config.bids_target_folder; 
cfg.participants                            = []; 

if isfield(config, 'session')
    cfg.ses                                     = config.session;
end
if isfield(config, 'run')
    cfg.run                                     = config.run;
end
if isfield(config, 'task')
    cfg.task                                    = config.task;
else
    cfg.task                                    = 'defaultTask'; 
end

% participant information
if exist('subjectInfo', 'var')
    
    allColumns      = subjectInfo.cols;
    
    % find the index of the subject nr column
    for iCol = 1:numel(allColumns)
        if strcmp(subjectInfo.cols(iCol), 'nr')
            nrColInd = iCol; 
        end
    end
    
    if ~exist('nrColInd','var')
        error('Participant info was provided without column "nr".')
    end
    
    % find the column that contains information from the given participant
    Pi = find([subjectInfo.data{:,nrColInd}] == config.subject); % find the matching participant number
    
    for iCol = 1:numel(allColumns)
            cfg.participants.(allColumns{iCol}) = subjectInfo.data{Pi, iCol};
    end
    
end

%% 
% load and assign streams (parts taken from xdf2fieldtrip)
%--------------------------------------------------------------------------
disp('Loading .xdf streams ...')
streams                  = load_xdf(cfg.dataset);

% initialize an array of booleans indicating whether the streams are continuous
iscontinuous = false(size(streams));

names = {};

% figure out which streams contain continuous/regular and discrete/irregular data
for i=1:numel(streams)
    
    names{i}           = streams{i}.info.name;
    
    % if the nominal srate is non-zero, the stream is considered continuous
    if ~strcmpi(streams{i}.info.nominal_srate, '0')
        
        num_samples  = numel(streams{i}.time_stamps);
        
        if num_samples > 20 % assume at least 20 samples in a continuous data stream

            iscontinuous(i) =  true;
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
            
        end
        
    else
        try
            num_samples  = numel(streams{i}.time_stamps);
            t_begin      = streams{i}.time_stamps(1);
            t_end        = streams{i}.time_stamps(end);
            duration     = t_end - t_begin;
            
            % if sampling rate is higher than 20 Hz,
            % the stream is considered continuous
            if (num_samples - 1) / duration >= 20 && duration > 1
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

if importEEG
    eegStreamName = config.eeg.stream_name; 
    xdfeeg        = streams(contains(lower(names),lower(eegStreamName)) & iscontinuous);
    
    if isempty(xdfeeg)
        error('No eeg streams found - check whether stream_name match the names of streams in .xdf')
    elseif numel(xdfeeg) > 1
        error('Multiple eeg streams found - usage not supported')
    end
end

if importMotion
    for Si = 1:numel(config.motion.streams) 
        motionStreamNames{Si}   = config.motion.streams{Si}.stream_name;
    end
    xdfmotion   = streams(contains(lower(names),lower(motionStreamNames)) & iscontinuous);
    
    if isempty(xdfmotion)
        error('Configuration field motion specified but no streams found - check whether stream_name match the names of streams in .xdf')
    end

end

if importPhys
    for Si = 1:numel(config.phys.streams)
        physioStreamNames{Si}   = config.phys.streams{Si}.stream_name;
    end
    xdfphysio   = streams(contains(lower(names),lower(physioStreamNames)) & iscontinuous);
    
    if isempty(xdfphysio)
        error('Configuration field physio specified but no streams found - check whether stream_name match the names of streams in .xdf')
    end
end

xdfmarkers  = streams(~iscontinuous);

%%
if importEEG % This loop is always executed in current version
    
    %----------------------------------------------------------------------
    %                   Convert EEG Data to BIDS
    %----------------------------------------------------------------------
    
    % construct fieldtrip data
    eeg        = stream2ft(xdfeeg{1});
    
    % save eeg start time
    eegStartTime                = eeg.time{1}(1);
    
    % eeg metadata construction
    %----------------------------------------------------------------------
    eegcfg                              = cfg;
    eegcfg.datatype                     = 'eeg';
    eegcfg.method                       = 'convert';
    
    % default coordinate system files 
    if isfield(config.eeg, 'chanloc')
        if ~isempty(config.eeg.chanloc)
            eegcfg.coordsystem.EEGCoordinateSystem      = 'n/a';
            eegcfg.coordsystem.EEGCoordinateUnits       = 'mm';
        end
    elseif isfield(config.eeg, 'elec_struct')
        if ~isempty(config.eeg.elec_struct)
            eegcfg.coordsystem.EEGCoordinateSystem      = 'n/a';
            eegcfg.coordsystem.EEGCoordinateUnits       = 'mm';
        end
    end
    
    % try to use information from preprocessing config
    if isfield(config.eeg, 'ref_channel')
        eegcfg.eeg.EEGReference                 = config.ref_channel; % field name comes from bemobil preprocessing pipeline
    end
    
    if isfield(config.eeg, 'linefreqs')
        if numel(config.linefreqs) == 1
            eegcfg.eeg.PowerLineFrequency           = config.linefreqs; % field name comes from bemobil preprocessing pipeline
        elseif numel(config.linefreqs) > 1
            eegcfg.eeg.PowerLineFrequency           = config.linefreqs(1);
            warning('Only the first value specified in config.eeg.linefreqs entered in eeg.json')
        end
    end
    
    % try to use metadata provided by the user - if provided, will overwrite values from config. 
    if exist('eegInfo','var')
        if isfield(eegInfo, 'eeg')
            eegcfg.eeg          = eegInfo.eeg;
        end
        if isfield(eegInfo, 'coordsystem')
            eegcfg.coordsystem  = eegInfo.coordsystem;
        end
    end
    
    % check if mandatory fields are specified and if not, fill with default values
    if isfield(eegcfg, 'eeg')
        [eegcfg.eeg] =  checkfield(eegcfg.eeg, 'EEGReference', 'REF', 'REF');
        [eegcfg.eeg] =  checkfield(eegcfg.eeg, 'PowerLineFrequency', 'n/a', 'n/a');
        [eegcfg.eeg] =  checkfield(eegcfg.eeg, 'SoftwareFilters', 'n/a', 'n/a');
    else
        eegcfg.eeg.EEGReference = 'REF';
        eegcfg.eeg.PowerLineFrequency = 'n/a';
        eegcfg.eeg.SoftwareFilters = 'n/a';

    end
    
    % read in the event stream (synched to the EEG stream)
    if ~isempty(xdfmarkers)
        
        if any(cellfun(@(x) ~isempty(x.time_series), xdfmarkers))
            
            events                  = stream2events(xdfmarkers, xdfeeg{1}.time_stamps);
            eventsFound             = 1;
            
            % event parser script
            if isfield(config, 'bids_parsemarkers_custom')
                if isempty(config.bids_parsemarkers_custom)
                    [events, eventsJSON] = bemobil_bids_parsemarkers(events);
                else
                    [events, eventsJSON] = feval(config.bids_parsemarkers_custom, events);
                end
            else 
                [events, eventsJSON] = bemobil_bids_parsemarkers(events);
            end
            
            eegcfg.events = events;
            
        end
    end
    
    if isfield(config.eeg, 'elec_struct')
        eegcfg.elec                         = config.eeg.elec_struct;
    elseif isfield(config.eeg, 'chanloc')
        if isfield(config.eeg, 'chanloc_newname')
            elec = ft_read_sens(config.eeg.chanloc);
            elec.label = config.eeg.chanloc_newname; 
            eegcfg.elec = elec; 
        else
            eegcfg.elec                         = config.eeg.chanloc;
        end
    end
    
    % acquisition time processing 
    eegcfg.acq_time = datestr(datenum(config.acquisition_time),'yyyy-mm-ddTHH:MM:SS.FFF'); % microseconds are rounded
    
    % write eeg files in bids format
    data2bids(eegcfg, eeg);
    
end

%%
if importMotion
    
    %----------------------------------------------------------------------
    %                   Convert Motion Data to BIDS
    %----------------------------------------------------------------------

    % construct motion metadata applying to all tracking systems
    %----------------------------------------------------------------------
    motioncfg                                         = cfg;                % copy general fields
    
    % data type
    motioncfg.datatype                                = 'motion';
    
    % coordinate system
    motioncfg.coordsystem.MotionCoordinateSystem      = motionInfo.coordsystem.MotionCoordinateSystem;
    motioncfg.coordsystem.MotionRotationRule          = motionInfo.coordsystem.MotionRotationRule;
    motioncfg.coordsystem.MotionRotationOrder         = motionInfo.coordsystem.MotionRotationOrder;
    motioncfg.TrackingSystemCount                     = numel(trackSysInData);


    % construct fieldtrip data
    ftmotion = {};
    for iM = 1:numel(xdfmotion)
        ftmotion{iM} = stream2ft(xdfmotion{iM});
    end
    
    % iterate over tracking systems
    %----------------------------------------------------------------------
    % initialize variables to concatenate over
    motioncfg.channels.name                 = {};
    motioncfg.channels.tracking_system      = {};
    motioncfg.channels.tracked_point        = {};
    motioncfg.channels.component            = {};
    motioncfg.channels.placement            = {};
    MotionChannelCount = 0;
    
    for tsi = 1:numel(trackSysInData)
        
        motionStreamNames   = kv_trsys_to_st(trackSysInData{tsi});
        trackedPointNames   = kv_trsys_to_trp(trackSysInData{tsi}); 
        anatomicalNames     = kv_trsys_to_anat(trackSysInData{tsi}); 
        
        % check if the names are nested in a cell
        newCell = {};
        for tpi = 1:numel(trackedPointNames)
            if iscell(trackedPointNames(tpi))
                newCell = [newCell trackedPointNames{tpi}];
            else
                newCell{end + 1} = trackedPointNames(tpi);
            end
        end
        trackedPointNames = newCell; 
        
        % check if the names are nested in a cell
        newCell = {};
        for tpi = 1:numel(anatomicalNames)
            if iscell(anatomicalNames(tpi))
                newCell = [newCell anatomicalNames{tpi}];
            else
                newCell{end + 1} = anatomicalNames(tpi);
            end
        end
        
        anatomicalNames = newCell;
        
        if isfield(cfg, 'ses')
            si = cfg.ses;
        else
            si = 1;
        end
        if isfield(cfg, 'run')
            ri = cfg.run;
        else
            ri = 1;
        end
        
        streamInds = [];
        for Fi = 1:numel(ftmotion)
           if contains(lower(ftmotion{Fi}.hdr.orig.name),lower(motionStreamNames))
               streamInds(end+1) = Fi; 
           end
        end
        
        % if needed, execute a custom function for any alteration to the data to address dataset specific issues
        % (quat2eul conversion, unwrapping of angles, resampling, wrapping back to [pi, -pi], and concatenating for instance)
        motion = feval(motionCustom, ftmotion(streamInds), trackedPointNames, config.subject, si, ri);
      
        % channel metadata 
        %------------------------------------------------------------------
        rb_streams = motionStreamNames;
        rb_names = trackedPointNames ;
        rb_anat = anatomicalNames;
        
        for ci  = 1:motion.hdr.nChans
            
            motionChanType          = motion.hdr.chantype{ci};
            
            if isfield(config.motion, motionChanType)
                motion.hdr.chanunit{ci} = config.motion.(motionChanType).unit; 
                disp(['Using custom unit ' config.motion.(motionChanType).unit ' for type ' motionChanType])
            else
                motion.hdr.chanunit{ci} = motion_type.(motionChanType).unit;
            end
            
            splitlabel                                      = regexp(motion.hdr.label{ci}, '_', 'split');
            motioncfg.channels.name{end+1}                  = motion.hdr.label{ci};
            motioncfg.channels.tracking_system{end+1}       = trackSysInData{tsi};
            
            % assign object names and anatomical positions
            for iN = 1:numel(rb_names)
                if contains(lower(motion.hdr.label{ci}),lower(rb_names{iN}))
                    motioncfg.channels.tracked_point{end+1}        = rb_names{iN};
                    motioncfg.channels.placement{end+1}            = rb_anat{iN};
                end
            end
            
            motioncfg.channels.component{end+1}    = splitlabel{end};
            
        end
        
        % tracking system-specific information 
        %------------------------------------------------------------------
        motioncfg.tracksys                                = trackSysInData{tsi};
        
        % sampling frequency
        if isfield(motion, 'fsample')
            effectiveSRate = motion.fsample;
        else
            effectiveSRate = motion.hdr.Fs;
        end
        
        % copy motion metadata fields
        motioncfg.motion = motionInfo.motion;
        
        % start time
        motionStartTime                 = motion.time{1}(1);
        motionTimeShift                 = motionStartTime - eegStartTime;
        
        % shift acq_time to store relative offset to eeg data
        acq_time = datenum(config.acquisition_time) + (motionTimeShift/(24*60*60));
        motioncfg.scans.acq_time = datestr(acq_time,'yyyy-mm-ddTHH:MM:SS.FFF'); % milisecond precision
  
        % effective sampling rate 
        motioncfg.motion.TrackingSystems.(trackSysInData{tsi}).SamplingFrequencyEffective = effectiveSRate; 
        
        % RecordingDuration
        motioncfg.motion.TrackingSystems.(trackSysInData{tsi}).RecordingDuration = (motion.hdr.nSamples*motion.hdr.nTrials)/effectiveSRate;

        % add the number of tracking points to tracked point count
        motioncfg.motion.TrackingSystems.(trackSysInData{tsi}).TrackedPointsCount = sum(numel(trackedPointNames)); % add entries which contain rb_name for corresponding tracking system
        
        % add the number of channels to MotionChannelCount
        motioncfg.motion.MotionChannelCount = MotionChannelCount + motion.hdr.nChans;
        
        % write motion files in bids format
        data2bids(motioncfg, motion);
    end

end

%%
if importPhys
    
    %----------------------------------------------------------------------
    %               Convert Generic Physio Data to BIDS
    %----------------------------------------------------------------------
    
    ftphysio = {};
    
    % construct fieldtrip data
    for iP = 1:numel(xdfphysio)
        ftphysio{iP} = stream2ft(xdfphysio{iP});
    end
    
    for Pi = 1:numel(PhyStreamsInData)
        
        % resample data to match the stream of highest srate (no custom processing supported for physio data yet)
        physio = feval(physioCustom, ftphysio, physioStreamNames, participantNr, si, di);
        
        % construct physio metadata
        physiocfg               = cfg;                                           % copy general fields
        physiocfg.datatype      = 'physio';
        
        %------------------------------------------------------------------
        if ~exist('physioInfo', 'var')
            
            % default values for physio specific fields in json
            physioInfo.physio.Manufacturer                     = 'Undefined';
            physioInfo.physio.ManufacturersModelName           = 'Undefined';
            physioInfo.physio.RecordingType                    = 'continuous';
            
        end
        
        % physio specific fields in json
        physiocfg.physio                                  = physioInfo.physio;
        
        % start time
        physioStartTime                 = physio.time{1}(1);
        physioTimeShift                 = physioStartTime - eegStartTime;
        
        % shift acq_time to store relative offset to eeg data
        acq_time = datenum(config.acquisition_time) + (physioTimeShift/(24*60*60));
        physiocfg.scans.acq_time = datestr(acq_time,'yyyy-mm-ddTHH:MM:SS.FFF');
        
        % start time
        physiocfg.physio.StartTime                        = physioStartTime - eegStartTime;
        
        % write motion files in bids format
        data2bids(physiocfg, physio);
        
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
% This function checks for existence of required fields 
% for optional fields that have default values, assign them 

% throw an error if a required field has not been specified
if strcmp(defaultValue, 'required')
    if ~isfield(oldconfig, fieldName)
        error(['Required config field ' fieldName ' not specified'])
    end
end

% assign default value to an optional field 
newconfig   = oldconfig;

if ~isfield(oldconfig, fieldName)
    newconfig.(fieldName) = defaultValue;
    warning(['Config field ' fieldName ' not specified- using default value: ' defaultValueText])
end


end

%--------------------------------------------------------------------------
function [ftdata] = stream2ft(xdfstream)

% construct header
hdr.Fs                  = xdfstream.info.effective_srate;
hdr.nFs                 = str2double(xdfstream.info.nominal_srate);
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

        try 
            hdr.chantype{j} = xdfstream.info.desc.channels.channel{j}.type;
        catch
            disp([hdr.label{j} ' missing type'])
        end
        
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
