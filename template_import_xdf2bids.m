
% initialize EEGLAB 
if ~exist('ALLCOM','var')
	eeglab;
end

% initialize fieldtrip without adding alternative files to path 
global ft_default
ft_default.toolbox.signal = 'matlab';  % can be 'compat' or 'matlab'
ft_default.toolbox.stats  = 'matlab';
ft_default.toolbox.image  = 'matlab';
ft_defaults % this sets up the FieldTrip path

%% [OPTIONAL] check the .xdf data to explore the structure
ftPath      = fileparts(which('ft_defaults'));
addpath(fullfile(ftPath, 'external','xdf')); 
xdfPath     = '...yourpath\yourfile.xdf'; % enter full path to your .xdf file 

% load .xdf data to check what is in there
streams         = load_xdf(xdfPath);
streamnames     = cellfun(@(x) x.info.name, streams, 'UniformOutput', 0)' % will display names of streams contained in .xdf

% display names of all channels in the .xdf data
for Si = 1:numel(streamnames)
    if isfield( streams{Si}.info.desc, 'channels')
        channelnames    = cellfun(@(x) x.label, streams{Si}.info.desc.channels.channel, 'UniformOutput', 0)'
    end
end


%% [OPTIONAL] enter metadata about the data set, data modalities, and participants

% general metadata shared across all modalities
% will be saved in BIDS-folder/data_description.json
% see "https://bids-specification.readthedocs.io/en/stable/03-modality-agnostic-files.html#:~:text=LICENSE-,dataset_description.json,-The%20file%20dataset_description"
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
generalInfo = [];

% required for dataset_description.json
generalInfo.dataset_description.Name                = 'name of your data set';
generalInfo.dataset_description.BIDSVersion         = 'version of BIDS-specification you are following';

% optional for dataset_description.json
generalInfo.dataset_description.License             = 'licence type';
generalInfo.dataset_description.Authors             = {"author 1", "author 2", "author 3"};
generalInfo.dataset_description.Acknowledgements    = 'acknowledgement text';
generalInfo.dataset_description.Funding             = {"funding source 1", "funding source 2"};
generalInfo.dataset_description.ReferencesAndLinks  = {"reference", "link to article"};
generalInfo.dataset_description.DatasetDOI          = 'DOI of your data set';

% general information shared across modality specific json files 
generalInfo.InstitutionName                         = 'name of your institute';
generalInfo.InstitutionalDepartmentName             = 'name of your department';
generalInfo.InstitutionAddress                      = 'address of your institute';
generalInfo.TaskDescription                         = 'text describing your task';
 

% information about the eeg recording system 
% will be saved in BIDS-folder/sub-XX/[ses-XX]/eeg/*_eeg.json and *_coordsystem.json
% see "https://bids-specification.readthedocs.io/en/stable/04-modality-specific-files/03-electroencephalography.html#:~:text=MAY%20be%20specified.-,Sidecar%20JSON%20(*_eeg.json),-Generic%20fields%20MUST"
% and "https://bids-specification.readthedocs.io/en/stable/04-modality-specific-files/03-electroencephalography.html#:~:text=after%20the%20recording.-,Coordinate%20System%20JSON,-(*_coordsystem.json"
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
eegInfo     = [];
eegInfo.coordsystem.EEGCoordinateSystem     = 'enter the name of your coordinate system'; % only needed when you share eloc
eegInfo.coordsystem.EEGCoordinateUnits      = 'enter the unit of your coordinate system'; % only needed when you share eloc
eegInfo.coordsystem.EEGCoordinateSystemDescription = 'enter description of your coordinate system'; % only needed when you share eloc
eegInfo.eeg.SamplingFrequency               = 1000; % nominal sampling frequency  
                                                  
% information about the motion recording system 
% will be saved in BIDS-folder/sub-XX/[ses-XX]/motion/*_motion.json 
% see "https://docs.google.com/document/d/1iaaLKgWjK5pcISD1MVxHKexB3PZWfE2aAC5HF_pCZWo/edit?usp=sharing"
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
motionInfo  = []; 

tracking_systems = {'System1', 'System2', 'System3'}; % enter the names of your tracking systems 

% motion specific fields in json
motionInfo.motion = [];

% system 1 information
motionInfo.motion.TrackingSystems(1).TrackingSystemName               = tracking_systems{1};
motionInfo.motion.TrackingSystems(1).Manufacturer                     = 'HTC';
motionInfo.motion.TrackingSystems(1).ManufacturersModelName           = 'Vive Pro';
motionInfo.motion.TrackingSystems(1).SamplingFrequency                = 90; %  If no nominal Fs exists, n/a entry returns 'n/a'. If it exists, n/a entry returns nominal Fs from motion stream.
motionInfo.motion.TrackingSystems(1).DeviceSerialNumber               = 'n/a'; 
motionInfo.motion.TrackingSystems(1).SoftwareVersions                 = 'n/a'; 
motionInfo.motion.TrackingSystems(1).SpatialAxes                      = 'FRU'; 
motionInfo.motion.TrackingSystems(1).RotationRule                     = 'left-hand'; 
motionInfo.motion.TrackingSystems(1).RotationOrder                    = 'ZXY'; 

% system 2 information
motionInfo.motion.TrackingSystems(2).TrackingSystemName               = tracking_systems{2};
motionInfo.motion.TrackingSystems(2).Manufacturer                     = 'Impuls X2';
motionInfo.motion.TrackingSystems(2).ManufacturersModelName           = 'PhaseSpace';
motionInfo.motion.TrackingSystems(2).SamplingFrequency                = 90;
motionInfo.motion.TrackingSystems(2).DeviceSerialNumber               = 'n/a'; 
motionInfo.motion.TrackingSystems(2).SoftwareVersions                 = 'n/a';
motionInfo.motion.TrackingSystems(2).SpatialAxes                       = 'FRU'; 
motionInfo.motion.TrackingSystems(2).RotationRule                     = 'left-hand'; 
motionInfo.motion.TrackingSystems(2).RotationOrder                    = 'ZXY'; 


% system 3 information
motionInfo.motion.TrackingSystems(3).TrackingSystemName               = tracking_systems{3};
motionInfo.motion.TrackingSystems(3).Manufacturer                     = 'Virtual System Manufacturer';
motionInfo.motion.TrackingSystems(3).ManufacturersModelName           = 'Virtual System Manufacturer Model';
motionInfo.motion.TrackingSystems(3).SamplingFrequency                = 60;
motionInfo.motion.TrackingSystems(3).DeviceSerialNumber               = 'n/a'; 
motionInfo.motion.TrackingSystems(3).SoftwareVersions                 = 'n/a'; 
motionInfo.motion.TrackingSystems(3).SpatialAxes                      = 'FRU'; 
motionInfo.motion.TrackingSystems(3).RotationRule                     = 'left-hand'; 
motionInfo.motion.TrackingSystems(3).RotationOrder                    = 'ZXY'; 


% participant information 
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
% here describe the fields in the participant file
% see "https://bids-specification.readthedocs.io/en/stable/03-modality-agnostic-files.html#participants-file:~:text=UTF%2D8%20encoding.-,Participants%20file,-Template%3A"
% for numerical values  : 
%       subjectData.fields.[insert your field name here].Description    = 'describe what the field contains';
%       subjectData.fields.[insert your field name here].Unit           = 'write the unit of the quantity';
% for values with discrete levels :
%       subjectData.fields.[insert your field name here].Description    = 'describe what the field contains';
%       subjectData.fields.[insert your field name here].Levels.[insert the name of the first level] = 'describe what the level means';
%       subjectData.fields.[insert your field name here].Levels.[insert the name of the Nth level]   = 'describe what the level means';
%--------------------------------------------------------------------------
subjectInfo.fields.nr.Description       = 'numerical ID of the participant'; 
subjectInfo.fields.age.Description      = 'age of the participant'; 
subjectInfo.fields.age.Unit             = 'years'; 
subjectInfo.fields.sex.Description      = 'sex of the participant'; 
subjectInfo.fields.sex.Levels.M         = 'male'; 
subjectInfo.fields.sex.Levels.F         = 'female'; 
subjectInfo.fields.handedness.Description    = 'handedness of the participant';
subjectInfo.fields.handedness.Levels.R       = 'right-handed';
subjectInfo.fields.handedness.Levels.L       = 'left-handed';

% names of the columns - 'nr' column is just the numerical IDs of subjects
%                         do not change the name of this column
subjectInfo.cols = {'nr',   'age',  'sex',  'handedness'};
subjectInfo.data = {1,     30,     'F',     'R' ; ...
                    2,     22,     'M',     'R'; ...
                    3,     23,     'F',     'R'; ...
                    4,     34,     'M',     'R'; ...
                    5,     25,     'F',     'R'; ...
                    6,     21,     'F',     'R' ; ...
                    7,     28,     'M',     'R'; ...
                    8,     28,     'M',     'R'; ...
                    9,     24,     'F',     'R'; ...
                    10,    25,     'F',     'L'; ...
                    11,    30,     'F',     'R'; ...
                    12,    22,     'M',     'R'; ...
                    13,    23,     'F',     'R'; ...
                    14,    34,     'M',     'R'; ...
                    15,    25,     'F',     'R'; ...
                    16,    21,     'F',     'R' ; ...
                    17,    28,     'M',     'R'; ...
                    18,    28,     'M',     'R'; ...
                    19,    24,     'F',     'R'; ...
                    20,    25,     'F',     'L';};
               

                
%%
% loop over participants
for subject = 1:20
    
    % loop over sessions 
    for session = 1:2
        
        config                        = [];                                 % reset for each loop ´
        config.bids_target_folder     = 'P:\Sein_Jeung\Project_BIDS\Example_datasets\SPOT_rotation\1_BIDS-data'; % required
        config.filename               = fullfile(['P:\Sein_Jeung\Project_BIDS\Example_datasets\SPOT_rotation\0_raw-data\' num2str(subject) '\test_' sessionNames{session} '.xdf']); % required
        config.eeg.chanloc            = fullfile(['P:\Sein_Jeung\Project_BIDS\Example_datasets\SPOT_rotation\0_raw-data\' num2str(subject) '\channel_locations.elc']); % optional
        
        config.task                   = 'Rotation';                         % optional
        config.subject                = subject;                            % required
        config.session                = sessionNames{session};              % optional
        config.overwrite              = 'on';
        
        config.eeg.stream_name        = 'BrainVision';                      % required
        % config.eeg. sampling frequency % override xdf nominal srate 
        
        %------------------------------------------------------------------
        
        if session == 1 % session body
            config.motion.streams{1}.xdfname            = 'headrigid';
            config.motion.streams{1}.bidsname           = tracking_systems{1};
            config.motion.streams{1}.tracked_points     = 'headRigid';
            config.motion.streams{1}.tracked_points_anat= 'head';
            config.motion.streams{1}.positions.channel_names = {'headRigid_Rigid_headRigid_X';  'headRigid_Rigid_headRigid_Y' ; 'headRigid_Rigid_headRigid_Z' }; 
            config.motion.streams{1}.quaternions.channel_names  = {'headRigid_Rigid_headRigid_quat_W';'headRigid_Rigid_headRigid_quat_Z';...
                                                                  'headRigid_Rigid_headRigid_quat_X';'headRigid_Rigid_headRigid_quat_Y'};
            
            if subject ~= 6 && subject ~= 13 && subject ~= 16 && subject ~= 19 && subject ~= 20 % subject ~= 6 && subject ~= 13 && subject ~= 19 && subject ~= 20 % subject 6 missing phasespace data and sub 13, 19, 20 has all zero stream
                config.motion.streams{2}.xdfname        = 'AllPhaseSpace';
                config.motion.streams{2}.bidsname    = tracking_systems{2};
                config.motion.streams{2}.tracked_points     = {'Rigid1', 'Rigid2', 'Rigid3', 'Rigid4'};
                config.motion.streams{2}.positions.channel_names = {'vizardAllPhasespaceLog_Rigid1_X', 'vizardAllPhasespaceLog_Rigid2_X', 'vizardAllPhasespaceLog_Rigid3_X', 'vizardAllPhasespaceLog_Rigid4_X';...
                                                                     'vizardAllPhasespaceLog_Rigid1_Y', 'vizardAllPhasespaceLog_Rigid2_Y', 'vizardAllPhasespaceLog_Rigid3_Y', 'vizardAllPhasespaceLog_Rigid4_Y';...
                                                                     'vizardAllPhasespaceLog_Rigid1_Z', 'vizardAllPhasespaceLog_Rigid2_Z', 'vizardAllPhasespaceLog_Rigid3_Z', 'vizardAllPhasespaceLog_Rigid4_Z'};
                config.motion.streams{2}.quaternions.channel_names = {'vizardAllPhasespaceLog_Rigid1_A', 'vizardAllPhasespaceLog_Rigid2_A', 'vizardAllPhasespaceLog_Rigid3_A', 'vizardAllPhasespaceLog_Rigid4_A';...
                                                                     'vizardAllPhasespaceLog_Rigid1_B', 'vizardAllPhasespaceLog_Rigid2_B', 'vizardAllPhasespaceLog_Rigid3_B', 'vizardAllPhasespaceLog_Rigid4_B';...
                                                                     'vizardAllPhasespaceLog_Rigid1_C', 'vizardAllPhasespaceLog_Rigid2_C', 'vizardAllPhasespaceLog_Rigid3_C', 'vizardAllPhasespaceLog_Rigid4_C'; ...
                                                                     'vizardAllPhasespaceLog_Rigid1_D', 'vizardAllPhasespaceLog_Rigid2_D', 'vizardAllPhasespaceLog_Rigid3_D', 'vizardAllPhasespaceLog_Rigid4_D'};
                                                              
            end
        else
            config.motion.streams{1}.xdfname                    = 'head';
            config.motion.streams{1}.bidsname                   = tracking_systems{3};
            config.motion.streams{1}.tracked_points             = 'headRigid';
            config.motion.streams{1}.positions.channel_names    = {'headRigid_Rigid_headRigid_X';  'headRigid_Rigid_headRigid_Y' ; 'headRigid_Rigid_headRigid_Z' }; 
            config.motion.streams{1}.quaternions.channel_names  = {'headRigid_Rigid_headRigid_quat_W';'headRigid_Rigid_headRigid_quat_Z';...
                                                                    'headRigid_Rigid_headRigid_quat_X';'headRigid_Rigid_headRigid_quat_Y'};
            config.motion.POS.unit                              = 'vm';
        end
        
        bemobil_xdf2bids(config, ...
            'general_metadata', generalInfo,...
            'participant_metadata', subjectInfo,...
            'motion_metadata', motionInfo, ...
            'eeg_metadata', eegInfo);
    end
    
    % configuration for bemobil bids2set
    %----------------------------------------------------------------------
    config.set_folder             = studyFolder;
    config.session_names            = sessionNames;
    config.raw_EEGLAB_data_folder   = '2_raw-EEGLAB-testing';
    
    % match labels in electrodes.tsv and channels.tsv
    matchlocs = {};
    letters = {'g', 'y', 'r', 'w', 'n'};
    for Li = 1:numel(letters)
        letter = letters{Li}; 
        for Ni = 1:32
            matchlocs{Ni + (Li-1)*32,1} = [letter num2str(Ni)]; % channel name in electrodes.tsv
            matchlocs{Ni + (Li-1)*32,2} = ['BrainVision RDA_' upper(letter) num2str(Ni, '%02.f')]; % channel name in channels.tsv 
        end
    end
    
    [matchlocs{157:159,1}] = deal(''); 
    config.match_electrodes_channels = matchlocs; 
    config.other_data_types = {'motion'}; 
    config.use_nominal_srate = {'phasespace'};
%     bemobil_bids2set(config);

end