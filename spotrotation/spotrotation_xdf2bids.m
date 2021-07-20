% example script for BeMoBIL BIDS tools xdf2bids

% setup paths
addpath('C:\Users\sgrot\Documents\Uni\01_Master\6. Semester\00_MA\Themenfindung\BIDs\tools\eeglab_current\eeglab2021.0');
addpath(genpath('C:\Users\sgrot\Documents\Uni\01_Master\6. Semester\00_MA\Themenfindung\BIDs\code\bemobil-pipeline'));
addpath('C:\Users\sgrot\Documents\Uni\01_Master\6. Semester\00_MA\Themenfindung\BIDs\tools\fieldtrip');
bemobil_config.study_folder             = 'C:\Users\sgrot\Documents\Uni\01_Master\6. Semester\00_MA\Themenfindung\BIDs\data\';
eeglab;
ft_defaults
ftPath = fileparts(which('ft_defaults'));
addpath(fullfile(ftPath, 'external','xdf')); 

% check 
which bemobil_xdf2bids 
which load_xdf

% streamNr = 3; 
% % eegStreamNr = 5; 
% 
% % load .xdf data to check what is in there
% % streams = load_xdf(fullfile(bemobil_config.study_folder, '\spotrotation\vp-6\vp-6_control_body.xdf'));
% streams = load_xdf(fullfile(bemobil_config.study_folder, '\spotrotation\vp-7\vp-7_control_joy.xdf'));
% streamnames     = cellfun(@(x) x.info.name, streams, 'UniformOutput', 0)'
% channelnames    = cellfun(@(x) x.label, streams{streamNr}.info.desc.channels.channel, 'UniformOutput', 0)'
% 
% % visualize position streams 
% figure; plot(streams{streamNr}.time_series(1:3,:)', 'LineWidth', 2)
% set(gca,'FontSize',15) 
% legend('X', 'Y', 'Z')
% title(['Stream ' streamnames{streamNr} ' position streams'], 'Interpreter','none')
% 
% % visualize orientation streams 
% figure; plot(streams{streamNr}.time_series(4:7,:)', 'LineWidth', 2)
% set(gca,'FontSize',15) 
% legend('A', 'B', 'C', 'D')
% title(['Stream ' streamnames{streamNr} ' quaternion streams'], 'Interpreter','none')



% configuration 
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
bemobil_config.filename_prefix          = 'vp-';
bemobil_config.source_data_folder       = 'spotrotation\';
bemobil_config.session_names            = {'body', 'joy'}; 
bemobil_config.other_data_types         = {'motion'};  
bemobil_config.resample_freq            = [];
bemobil_config.channel_locations_filename = [];

% these streams should be processed as rigid body streams containing 3 dof position and 3 dof orientation data (e.g. derivatives and filters applied)
bemobil_config.rigidbody_streams        = {'headRigid', 'headRigid'};
bemobil_config.rigidbody_names          =  {'Head', 'JoyStickRotation'}; 
bemobil_config.rigidbody_anat           = {'head', 'n/a'}; 

% motion data 
bemobil_config.bids_rb_in_sessions                 = [1, 0; 0, 1];

bemobil_config.bids_motion_position_units              = {'m', 'virtual meters'};          % if multisession, cell array of size 1 x session number
bemobil_config.bids_motion_orientation_units           = {'rad', 'radians'};               % if multisession, cell array of size 1 x session number
bemobil_config.bids_motion_velocity_units              = {'m/s', 'm/s'};                   % if multisession, cell array of size 1 x session number
bemobil_config.bids_motion_angularvelocity_units       = {'rad/s', 'rad/s'};               % if multisession, cell array of size 1 x session number
bemobil_config.bids_motion_acceleration_units          = {'m/s^2', 'm/s^2'};               % if multisession, cell array of size 1 x session number
bemobil_config.bids_motion_angularacceleration_units   = {'rad/s^2', 'rad/s^2'};           % if multisession, cell array of size 1 x session number
bemobil_config.bids_motion_mangeticfield_units         = {'T', 'T'};                       % if multisession, cell array of size 1 x session number
bemobil_config.bids_motion_jointangle_units            = {'rad', 'rad'};                   % if multisession, cell array of size 1 x session number

bemobil_config.bids_data_folder         = '1_BIDS-data\';
bemobil_config.bids_eeg_keyword          = 'BrainVision RDA';                  % marker streams also contain these strings. However, only the continuous stream is imported
bemobil_config.bids_tasklabel           = 'spotrotation';

% custom function names - customization recommended for data sets that have
%                         an 'unconventional' naming scheme for motion channels
bemobil_config.bids_motionconvert_custom    = 'spotrotation_motionconvert';
bemobil_config.bids_parsemarkers_custom     = [];

% general metadata shared across all modalities
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
generalInfo = [];

% root directory (where you want your bids data to be saved)
generalInfo.bidsroot                                = fullfile(bemobil_config.study_folder, bemobil_config.bids_data_folder); 

% required for dataset_description.json
generalInfo.dataset_description.Name                = 'Spot rotation with body and joystick';
generalInfo.dataset_description.BIDSVersion         = 'unofficial extension';

% optional for dataset_description.json
generalInfo.dataset_description.License             = 'n/a';
generalInfo.dataset_description.Authors             = 'Gramann, K., Hohlefeld, F.U., Gehrke, L., Klug, M';
generalInfo.dataset_description.Acknowledgements    = 'Acknowledgements here';
generalInfo.dataset_description.Funding             = 'n/a';
generalInfo.dataset_description.ReferencesAndLinks  = 'n/a';
generalInfo.dataset_description.DatasetDOI          = 'n/a';

% general information shared across modality specific json files 
generalInfo.InstitutionName                         = 'Technische Universitaet zu Berlin';
generalInfo.InstitutionalDepartmentName             = 'Biological Psychology and Neuroergonomics';
generalInfo.InstitutionAddress                      = 'Strasse des 17. Juni 135, 10623, Berlin, Germany';
generalInfo.TaskDescription                         = 'Participants equipped with VR HMD rotated either physically or using a joystick.';
generalInfo.task                                    = bemobil_config.bids_tasklabel;  

% acquisition dates
                                                         %rows and columns corresponding to dates below
                                              
                                                         %             |   eeg   |    motion
                                                         %-----|-------|---------|-------------
                                                         %sub6 | body  |         |
                                                         %sub6 | joy   |         |
                                                         %sub7 | body  |         |
                                                         %sub7 | joy   |         |
                                             
generalInfo.dateList                                = {'2021-07-09 12:15:00' '2021-07-09 12:15:00' ...
                                                       '2021-07-10 12:15:00' '2021-07-10 12:15:00' ...
                                                       '2021-07-11 12:15:00' '2021-07-11 12:15:00' ...
                                                       '2021-07-12 12:15:00' '2021-07-12 12:15:00'} ;
                           
% shifting factor
generalInfo.shiftFac                                = randi([-1000,1000]);


% information about the motion recording system 
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
% data type and acquisition label
motionInfo.acq                                     = [];

% motion specific fields in json
motionInfo.motion = [];
motionInfo.motion.RecordingType                       = 'continuous';
motionInfo.motion.TrackingSystems                     = []; 
motionInfo.motion.tracksys                            = {'OPTpos', 'VIRpos'}; % index corresponds to sessions number 

% system 1 information

motionInfo.motion.TrackingSystems.(motionInfo.motion.tracksys{1}).Manufacturer                     = 'HTC';
motionInfo.motion.TrackingSystems.(motionInfo.motion.tracksys{1}).ManufacturersModelName           = 'Vive Pro';
motionInfo.motion.TrackingSystems.(motionInfo.motion.tracksys{1}).SamplingFrequencyNominal         = 'n/a'; %  If no nominal Fs exists, n/a entry returns 'n/a'. If it exists, n/a entry returns nominal Fs from motion stream.
motionInfo.motion.TrackingSystems.(motionInfo.motion.tracksys{1}).SamplingFrequencyEffective       = [];
motionInfo.motion.TrackingSystems.(motionInfo.motion.tracksys{1}).SoftwareFilters                  = 'n/a';

% system 2 information
motionInfo.motion.TrackingSystems.(motionInfo.motion.tracksys{2}).Manufacturer                     = 'Virtual System Manufacturer';
motionInfo.motion.TrackingSystems.(motionInfo.motion.tracksys{2}).ManufacturersModelName           = 'Virtual System Manufacturer Model';
motionInfo.motion.TrackingSystems.(motionInfo.motion.tracksys{2}).SamplingFrequencyNominal         = 'n/a';
motionInfo.motion.TrackingSystems.(motionInfo.motion.tracksys{2}).SamplingFrequencyEffective       = [];
motionInfo.motion.TrackingSystems.(motionInfo.motion.tracksys{2}).SoftwareFilters                  = 'n/a';

% coordinate system
motionInfo.coordsystem.MotionCoordinateSystem      = 'RUF';
motionInfo.coordsystem.MotionRotationRule          = 'left-hand';
motionInfo.coordsystem.MotionRotationOrder         = 'ZXY';

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
subjectInfo.fields.handedness.Description    = 'handedness of the participant';
subjectInfo.fields.handedness.Levels.R       = 'right-handed';
subjectInfo.fields.handedness.Levels.L       = 'left-handed';

% names of the columns - 'nr' column is just the numerical IDs of subjects
%                         do not change the name of this column
subjectInfo.cols = {'nr',   'age',  'sex', 'handedness'};
subjectInfo.data = {6,     20,     'F',    'R' ; ...
                    7,     20,     'M',    'R' };
               

% numerical IDs 
%--------------------------------------------------------------------------
numericalIDs                            = [6,7]; 

bemobil_xdf2bids(bemobil_config, numericalIDs, 'general_metadata', generalInfo, 'motion_metadata', motionInfo, 'participant_metadata', subjectInfo)