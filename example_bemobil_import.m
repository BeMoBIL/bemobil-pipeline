% This example was used to import data set used in MoBI workshop, SanDiego 2022 
% All paths are examples and to be adjusted to match where you have your
% toolboxes and data on your PC
% TODO : where to find the data
%--------------------------------------------------------------------------
% author : sein.jeung@campus.tu-berlin.de

addpath('M:\8_Conferences\MoBI_Workshop_2022\public\toolboxes\eeglab2022.0')
addpath('M:\8_Conferences\MoBI_Workshop_2022\public\toolboxes\fieldtrip-master') 

% start eeglab to add bemobil pipeline to matlab path
eeglab;

% add load_xdf for inspecting data 
ftPath = fileparts(which('ft_defaults'));
addpath(fullfile(ftPath, 'external','xdf')); 

% data folder 
dataFolder                  = 'P:\Sein_Jeung\Project_BIDS\Workshop_2022'; 


%% [optional] inspect.xdf data
% file name 
xdfFilePath                 = fullfile(dataFolder, 'source-data\vp_24\vp_24_walk.xdf');

% load xdf streams 
xdfStreams                  = load_xdf(xdfFilePath); 
nStreams                    = numel(xdfStreams); 

% list all channel names per stream
for iStream = 1:nStreams
   
    xdfStreams{iStream}.info.name
    
    if isfield(xdfStreams{iStream}.info.desc, 'channels')
        cellfun(@(x) x.label, xdfStreams{iStream}.info.desc.channels.channel, 'UniformOutput', false)'
    end
    
end

% find out the right permutation for quaternion conversion 
bemobil_inspect_quaternion_convert(xdfStreams{1}, [4:7])


%% specify metadata
% general metadata shared across all modalities
%--------------------------------------------------------------------------
generalInfo = [];

% root directory (where you want your bids data to be saved)
generalInfo.bidsroot                                = fullfile(dataFolder, config.bids_target_folder);

% required for dataset_description.json
generalInfo.dataset_description.Name                = 'Walking task in the young and old';
generalInfo.dataset_description.BIDSVersion         = 'unofficial extension';

% optional for dataset_description.json
generalInfo.dataset_description.License             = 'CC BY 4.0';
generalInfo.dataset_description.Authors             = 'Janna Protzak and Klaus Gramann';
generalInfo.dataset_description.Acknowledgements    = 'n/a';
generalInfo.dataset_description.Funding             = 'This Study was realized by funding from the Federal Ministry of Education and Research (BMBF)';
generalInfo.dataset_description.ReferencesAndLinks  = 'n/a';
generalInfo.dataset_description.DatasetDOI          = 'n/a';

% general information shared across modality specific json files
generalInfo.InstitutionName                         = 'Technische Universitaet Berlin';
generalInfo.InstitutionalDepartmentName             = 'Junior research group FANS (Pedestrian Assistance System for Older Road Users), Department of Psychology and Ergonomics';
generalInfo.InstitutionAddress                      = 'Marchstr. 23, 10587, Berlin, Germany';
generalInfo.TaskDescription                         = 'Younger and older adults performed a visual discrimination task (button presses to peripheral presented LED flashes) during walking. Visual targets were either presented with or without preceding vibro-tactile cues';
generalInfo.task                                    = 'VisualDiscrimination';

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
subjectInfo.fields.nr.Description       = 'numerical ID of the participant';
subjectInfo.fields.age.Description      = 'age of the participant';
subjectInfo.fields.age.Unit             = 'years';
subjectInfo.fields.sex.Description      = 'sex of the participant';
subjectInfo.fields.sex.Levels.M         = 'male';
subjectInfo.fields.sex.Levels.F         = 'female';
subjectInfo.fields.group.Description    = 'experiment group';
subjectInfo.fields.group.Levels.young   = 'younger participants under 35';
subjectInfo.fields.group.Levels.old     = 'older participants over 65';
subjectInfo.fields.handedness.Description    = 'handedness of the participant';
subjectInfo.fields.handedness.Levels.R       = 'right-handed';
subjectInfo.fields.handedness.Levels.L       = 'left-handed';

% names of the columns - 'nr' column is just the numerical IDs of subjects
%                         do not change the name of this column

subjectInfo.cols = {'nr',   'age',  'sex',  'group',    'handedness'};
subjectInfo.data = {...
    24,     30,     'F',    'young',    'R' ; ...
    64,     71,     'F',    'old',      'R' ; ...
    66,     67,     'M',    'old',      'R' ; ...
    76,     34,     'M',    'young',	'R' ; ...
    78,     33,     'M',    'young',	'R' ; ...
    };

% information about the eeg recording system 
%--------------------------------------------------------------------------
eegInfo                                             = [];
eegInfo.eeg.PowerLineFrequency                      = 50;
eegInfo.eeg.EEGReference                            = 'FCz';
eegInfo.eeg.SamplingFrequency                       = 1000; 
eegInfo.coordsystem.EEGCoordinateSystem             = 'Other';
eegInfo.coordsystem.EEGCoordinateUnits              = 'mm';
eegInfo.coordsystem.EEGCoordinateSystemDescription  = 'Measured with Xensor';
                                                  
% information about the motion recording system 
%--------------------------------------------------------------------------
motionInfo  = []; 
motionInfo.motion = [];

% iterate over all tracking systems
for TSi=1:9
    motionInfo.motion.TrackingSystems(TSi).TrackingSystemName          = bidsNames{TSi};
    motionInfo.motion.TrackingSystems(TSi).Manufacturer                = 'Impuls X2';
    motionInfo.motion.TrackingSystems(TSi).ManufacturersModelName      = 'PhaseSpace';
    motionInfo.motion.TrackingSystems(TSi).SamplingFrequency           = 480;
    motionInfo.motion.TrackingSystems(TSi).DeviceSerialNumber          = 'n.a.';
    motionInfo.motion.TrackingSystems(TSi).SoftwareVersions            = 'n.a.';
    motionInfo.motion.TrackingSystems(TSi).SpatialAxes                 = 'RUF';
    motionInfo.motion.TrackingSystems(TSi).RotationRule                = 'left-hand';
    motionInfo.motion.TrackingSystems(TSi).RotationOrder               = 'ZXY';
end

%% fill out config fields common to all .xdf files 
config                        = [];
config.bids_target_folder     = '1_BIDS-data';                              % required, string, bids target folder to be created
config.task                   = 'VisualDiscrimination';                     % optional, string, task name, default value 'defaultTask'

config.eeg.stream_name        = 'BrainVision';                              % required, string, a unique keyword in EEG stream to be searched for

% cell array of stream and tracking system names to be assigned to the streams
xdfNames                      = {'PhaseSpace_Rigid1',... % stream name (or unique keyword) in .xdf file
                                 'PhaseSpace_Rigid2',...
                                 'PhaseSpace_Rigid3',...
                                 'PhaseSpace_Rigid4',...
                                 'PhaseSpace_Rigid5',...
                                 'PhaseSpace_Rigid6',...
                                 'PhaseSpace_Rigid7',...
                                 'PhaseSpace_Rigid8',...
                                 'PhaseSpace_Rigid9'};
bidsNames                     = {'PhaseSpaceHead',... % how the tracking system (= stream) should be renamed in BIDS
                                 'PhaseSpaceLeftThigh',...
                                 'PhaseSpaceLeftLowerLeg',...
                                 'PhaseSpaceLeftAnkle',...
                                 'PhaseSpaceLeftForeFoot',...
                                 'PhaseSpaceRightThigh',...
                                 'PhaseSpaceRightLowerLeg',...
                                 'PhaseSpaceRightAnkle',...
                                 'PhaseSpaceRightForeFoot'};

% iterate over streams to configure each of them
for Si = 1:9
    config.motion.streams{Si}.xdfname        = xdfNames{Si}; % required, keyword in stream name, searched for in field "xdfdata{streamIndex}.info.name"
    config.motion.streams{Si}.bidsname       = bidsNames{Si}; % optional, name to be assgined in BIDS file name key-value pair as tracking system name
    config.motion.streams{Si}.tracked_points = {['Rigid' num2str(Si)]}; % reuired, keyword in channel names, indicating which object (tracked point) is included in the stream
    config.motion.streams{Si}.positions.channel_names =  {['PhaseSpace_Rigid' num2str(Si) '_Rigid' num2str(Si) '_X'];...
        ['PhaseSpace_Rigid' num2str(Si) '_Rigid' num2str(Si) '_Y'];...
        ['PhaseSpace_Rigid' num2str(Si) '_Rigid' num2str(Si) '_Z']};
    config.motion.streams{Si}.quaternions.channel_names = {['PhaseSpace_Rigid' num2str(Si) '_Rigid' num2str(Si) '_A'];...
        ['PhaseSpace_Rigid' num2str(Si) '_Rigid' num2str(Si) '_B'];...
        ['PhaseSpace_Rigid' num2str(Si) '_Rigid' num2str(Si) '_C'];...
        ['PhaseSpace_Rigid' num2str(Si) '_Rigid' num2str(Si) '_D']};
end


%% iterate over participants and sessions to convert .xdf to BIDS formatted data
subjects        = cell2mat(subjectInfo.data(:,1))'; 
sessionNames 	= {'walk', 'stand'}; 

for subject = subjects
     
    config.subject                = subject;                                         % required, subject numerical ID
       
    for iSes = 1:2
        
        config.filename               = fullfile(dataFolder, ['source-data\vp_' num2str(subject) '\vp_' num2str(subject) '_' sessionNames{iSes} '.xdf']);
        config.session                = sessionNames{iSes};                                     % optional, string, session name if there were multiple sessions

        bemobil_xdf2bids(config, ...
            'general_metadata', generalInfo,...
            'participant_metadata', subjectInfo,...
            'motion_metadata', motionInfo, ...
            'eeg_metadata', eegInfo);
    end
    
    %% convert bids to set
    %----------------------------------------------------------------------
    config                          = [];
    config.set_folder               = dataFolder;
    config.session_names            = sessionNames;
    config.raw_EEGLAB_data_folder   = '1_raw-EEGLAB';
    config.other_data_types         = {'motion'};
    
    bemobil_bids2set(config);
end



