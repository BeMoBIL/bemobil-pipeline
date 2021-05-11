% general metadata shared across all modalities
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
general_info = [];

% root directory (where you want your bids data to be saved)
general_info.bidsroot                                = fullfile(bemobil_config.study_folder, bemobil_config.bids_data_folder); 

% required for dataset_description.json
general_info.dataset_description.Name                = 'Walking task in the young and old';
general_info.dataset_description.BIDSVersion         = 'unofficial extension';

% optional for dataset_description.json
general_info.dataset_description.License             = 'CC BY 4.0';
general_info.dataset_description.Authors             = 'Janna Protzak and Klaus Gramann';
general_info.dataset_description.Acknowledgements    = 'n/a';
general_info.dataset_description.Funding             = 'This Study was realized by funding from the Federal Ministry of Education and Research (BMBF)';
general_info.dataset_description.ReferencesAndLinks  = 'n/a';
general_info.dataset_description.DatasetDOI          = 'n/a';

% general information shared across modality specific json files 
general_info.InstitutionName                         = 'Technische Universitaet Berlin';
general_info.InstitutionalDepartmentName             = 'Junior research group FANS (Pedestrian Assistance System for Older Road Users), Department of Psychology and Ergonomics';
general_info.InstitutionAddress                      = 'Marchstr. 23, 10587, Berlin, Germany';
general_info.TaskDescription                         = 'Younger and older adults performed a visual discrimination task (button presses to peripheral presented LED flashes) during walking. Visual targets were either presented with or without preceding vibro-tactile cues';
general_info.task                                    = bemobil_config.bids_task_label;  

% information about the EEG recording system 
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
% coordinate system
eeg_info.eeg.PowerLineFrequency                       = 50;
eeg_info.eeg.EEGReference                            = 'FCz'; 

% information about the motion recording system 
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
% data type and acquisition label
motion_info.acq                                     = 'PhS';

% motion specific fields in json
motion_info.motion.Manufacturer                     = 'PhaseSpace';
motion_info.motion.ManufacturersModelName           = 'ImpulseX2';
motion_info.motion.RecordingType                    = 'continuous';

% coordinate system
motion_info.coordsystem.MotionCoordinateSystem      = 'RUF';
motion_info.coordsystem.MotionRotationRule          = 'left-hand';
motion_info.coordsystem.MotionRotationOrder         = 'ZXY'; % this follows the result after quat to euler conversion

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
subject_info.fields.age.Description      = 'age of the participant'; 
subject_info.fields.age.Unit             = 'years'; 
subject_info.fields.sex.Description      = 'sex of the participant'; 
subject_info.fields.sex.Levels.M         = 'male'; 
subject_info.fields.sex.Levels.F         = 'female'; 
subject_info.fields.group.Description    = 'experiment group';
subject_info.fields.group.Levels.young   = 'younger participants under 35';
subject_info.fields.group.Levels.old     = 'older participants over 65';
subject_info.fields.handedness.Description    = 'handedness of the participant';
subject_info.fields.handedness.Levels.R       = 'right-handed';
subject_info.fields.handedness.Levels.L       = 'left-handed';

% names of the columns - 'nr' column is just the numerical IDs of subjects
%                         do not change the name of this column
subject_info.cols = {'nr',   'age',  'sex',  'group',    'handedness'};
subject_info.data = {64,     71,     'F',    'old',      'R' ; ...
                    66,     67,     'M',    'old',      'R' ; ...
                    76,     34,     'M',    'young',    'R' ; ...
                    78,     33,     'M',    'young',    'R' ...
                    };

