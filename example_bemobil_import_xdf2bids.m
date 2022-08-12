bidsFolder      = '...\bids-data'; % root directory (where you want your bids data to be saved)
taskName        = 'Task'; % name of your task  
sessionNames    = {'ses1', 'ses1'}; % session names 

% for motion data - stream names described in .xdf file 
xdfStreamNames     = { 'Stream1',...
                       'Stream2',...
                       'Stream3'};
% rename streams to save them as separate tracking systems in BIDS
% easily understandable names are recommended
bidsTracksysNames  = { 'HeadMountedDisplay',...
                       'ControllerL',...
                       'ControllerR'};
                   
% names of the tracked points in each stream
% the keyword has to be included in the channel labels in .xdf 
trackedPoints      = { 'head', ...
                       'lefthand', ...
                       'righthand'}; 
                   

% general metadata shared across all modalities
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
generalInfo = [];
generalInfo.bidsroot                                = bidsFolder;  

% required for dataset_description.json
generalInfo.dataset_description.Name                = 'Data set name';
generalInfo.dataset_description.BIDSVersion         = '1.2.0'; % if sharing motion data, use "unofficial extension"

% optional for dataset_description.json
generalInfo.dataset_description.License             = 'CC BY 4.0';
generalInfo.dataset_description.Authors             = {"Author1, A.", "Author2, A."};
generalInfo.dataset_description.Acknowledgements    = 'n/a';
generalInfo.dataset_description.Funding             = 'This Study was realized by funding from...';
generalInfo.dataset_description.ReferencesAndLinks  = {""};
generalInfo.dataset_description.DatasetDOI          = {"publication title", "doi"};

% general information shared across modality specific json files 
generalInfo.InstitutionName                         = 'Your Universitry';
generalInfo.InstitutionalDepartmentName             = 'Department of ...';
generalInfo.InstitutionAddress                      = 'Address, city, country';
generalInfo.TaskDescription                         = 'description of the task';
generalInfo.task                                    = taskName;  

% information about the EEG recording system
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
% coordinate system
eegInfo                                             = [];
eegInfo.eeg.PowerLineFrequency                      = 50;
eegInfo.eeg.EEGReference                            = 'FCz';
eegInfo.coordsystem.EEGCoordinateSystem             = 'Other';
eegInfo.coordsystem.EEGCoordinateUnits              = 'mm';
eegInfo.coordsystem.EEGCoordinateSystemDescription  = 'Measured with Xensor';
eegInfo.eeg.SamplingFrequency                       = 1000;

% information about the motion recording system 
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
motionInfo  = [];

% % motion specific fields in json
motionInfo.motion = [];
motionInfo.motion.RecordingType                     = 'continuous';

% iterate over streams (= tracking systems) and enter metadata
for idx_Mi=1:3
    
    motionInfo.motion.TrackingSystems(idx_Mi).TrackingSystemName          = bidsTracksysNames{idx_Mi};
    motionInfo.motion.TrackingSystems(idx_Mi).Manufacturer                = '...'; % manufacturer of the motion capture system
    motionInfo.motion.TrackingSystems(idx_Mi).ManufacturersModelName      = '...'; % model name of the tracking system
    motionInfo.motion.TrackingSystems(idx_Mi).SamplingFrequency           = 480; % nomnial sampling frequency designated by manufacturer
    motionInfo.motion.TrackingSystems(idx_Mi).SpatialAxes                 = 'RUF'; % XYZ spatial axes description 
    motionInfo.motion.TrackingSystems(idx_Mi).RotationRule                = 'left-hand'; % rotation rule 
    motionInfo.motion.TrackingSystems(idx_Mi).RotationOrder               = 'ZXY'; % rotation order (external) for euler angles
    
end

% participant information 
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
% here describe the fields in the participant file
% for numerical values  : 
%       subjectData.fields.[insert your field name here].Description    = 'describe what the field contains';
%       subjectData.fields.[insert your field name here].Unit           = 'write the unit of the quantity';
% for values with discrete levels :
%       subjectData.fields.[insert your field name here].Description    = 'describe what the field contains';
%       subjectData.fields.[insert your field name here].Levels.[insert the name of the first level] = 'describe what the level means';
%       subjectData.fields.[insert your field name here].Levels.[insert the name of the Nth level]   = 'describe what the level means';
%--------------------------------------------------------------------------
subjectInfo.fields.age.Description      = 'age of the participant'; 
subjectInfo.fields.age.Unit             = 'years'; 
subjectInfo.fields.sex.Description      = 'sex of the participant'; 
subjectInfo.fields.sex.Levels.M         = 'male'; 
subjectInfo.fields.sex.Levels.F         = 'female'; 
subjectInfo.fields.group.Description    = 'experiment group';
subjectInfo.fields.group.Levels.young   = 'younger participants under ...'; % custom column
subjectInfo.fields.group.Levels.old     = 'older participants over ...'; % custom column
subjectInfo.fields.handedness.Description    = 'handedness of the participant';
subjectInfo.fields.handedness.Levels.R       = 'right-handed';
subjectInfo.fields.handedness.Levels.L       = 'left-handed';

% names of the columns - 'nr' column is just the numerical IDs of subjects
%                         do not change the name of this column
subjectInfo.cols = {'nr',   'age',  'sex',  'group',    'handedness'};
subjectInfo.data = {1,     71,     'F',    'old',      'R' ; ...
                    2,     67,     'M',    'old',      'R' ; ...
                    3,     34,     'M',    'young',    'R' ; ...
                    4,     33,     'M',    'young',    'R' ...
                    };
                
                
%% General Setup

% enter all subjects to process here (you can split it up in more MATLAB instances if you have more CPU power and RAM)
subjects = [1 2 3 4];


%%  loop over xdf files 
for subject = subjects
    
    % loop over sessions
    for session = 1:2
        
        config                        = []; % reset for each loop
        config.bids_target_folder     = bidsFolder; % required
        config.task                   = taskName;                           % optional
        config.subject                = subject;                            % required
        config.session                = sessionNames;                       % optional
        config.filename               = fullfile(['...\subject-' num2str(subject) '.xdf']); % required, the name of the .xdf file
        
        % eeg configuration
        config.eeg.stream_name        = 'EEGsystemName';                    % required, name of the eeg stream as specified in .xdf
       
        % Tracking System Description:
        
        % if you do not know the correct way to transform your quaternion motion tracking values to euler angles, this
        % script might help:
        
        % help bemobil_inspect_quaternion_convert

        % Stream description        
        for idx_str = 1:3
            config.motion.streams{idx_str}.xdfname            = xdfStreamNames{idx_str}; % required, keyword in stream name, searched for in field "xdfdata{streamIndex}.info.name"
            config.motion.streams{idx_str}.bidsname           = bidsTracksysNames{idx_str}; % required, match with one of the values in "motion.tracksys{}.name"
            config.motion.streams{idx_str}.tracked_points     = trackedPoints{idx_str};
            
            % names of position and quaternion channels in each stream
            % for simplicity assume that all 3 streams have same style of channel names
            % (tracked point name + channel suffix, for instance, "head_POS_X")
            config.motion.streams{idx_str}.positions.channel_names =  {[trackedPoints{idx_str}, '_POS_X'];...
                                                                        [trackedPoints{idx_str}, '_POS_Y']; ...
                                                                        [trackedPoints{idx_str}, '_POS_Z']};
            config.motion.streams{idx_str}.quaternions.channel_names =  {[trackedPoints{idx_str}, '_QUAT_W'];...
                                                                        [trackedPoints{idx_str}, '_QUAT_X']; ...
                                                                        [trackedPoints{idx_str}, '_QUAT_Y']; ...
                                                                        [trackedPoints{idx_str}, '_QUAT_Z']};        
        end

        %% run xdf2bids for this xdf file 
        bemobil_xdf2bids(config, ...
            'general_metadata', generalInfo,...
            'participant_metadata', subjectInfo,...
            'motion_metadata', motionInfo, ...
            'eeg_metadata', eegInfo);
        
    end
    
end


