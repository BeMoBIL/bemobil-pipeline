% example script for BeMoBIL BIDS tools xdf2bids

% configuration 
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
bemobil_config.study_folder             = 'M:\8_Conferences\MoBI Workshop\data\';
bemobil_config.filename_prefix          = 'vp_';
bemobil_config.raw_data_folder          = '0_raw-data\';
bemobil_config.filenames                = {'_walk'}; 
bemobil_config.resample_freq            = 250;

% these streams should be processed as rigid body streams containing 3 dof position and 3 dof orientation data (e.g. derivatives and filters applied)
bemobil_config.rigidbody_streams        = {'PhaseSpace_Rigid1','PhaseSpace_Rigid2', 'PhaseSpace_Rigid3', ...
                                            'PhaseSpace_Rigid4','PhaseSpace_Rigid5', 'PhaseSpace_Rigid6', ...
                                            'PhaseSpace_Rigid7','PhaseSpace_Rigid8', 'PhaseSpace_Rigid9'};
bemobil_config.channel_locations_filename = [];                            

% bids fields 
bemobil_config.bids_data_folder         = '1_BIDS-data\';
bemobil_config.bids_rbsessions          = true(1,numel(bemobil_config.rigidbody_streams));
bemobil_config.bids_eegkeyword          = {'BrainVision'};                  % marker streams also contain these strings. However, only the continuous stream is imported
bemobil_config.bids_tasklabel           = 'walking';

% optional custom function names - customization highly recommeded
bemobil_config.bids_motioncfg_custom    = 'bids_motioncfg_mobiworkshop';
bemobil_config.bids_motionconvert_custom  = 'bids_motionconvert_mobiworkshop';
bemobil_config.bids_eegcfg_custom       = [];
bemobil_config.bids_parsemarkers_custom = 'bids_parsemarkers_mobiworkshop';


% general metadata shared across all modalities
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
generalinfo = [];

% root directory (where you want your bids data to be saved)
generalinfo.bidsroot                                = fullfile(bemobil_config.study_folder, bemobil_config.bids_data_folder); 

% required for dataset_description.json
generalinfo.dataset_description.Name                = 'Walking task in the young and old';
generalinfo.dataset_description.BIDSVersion         = 'unofficial extension';

% optional for dataset_description.json
generalinfo.dataset_description.License             = 'n/a';
generalinfo.dataset_description.Authors             = 'Janna Protzak and Klaus Gramann';
generalinfo.dataset_description.Acknowledgements    = 'Acknowledgements here';
generalinfo.dataset_description.Funding             = 'n/a';
generalinfo.dataset_description.ReferencesAndLinks  = 'n/a';
generalinfo.dataset_description.DatasetDOI          = 'n/a';

% general information shared across modality specific json files 
generalinfo.InstitutionName                         = 'Technische Universitaet zu Berlin';
generalinfo.InstitutionalDepartmentName             = 'Biological Psychology and Neuroergonomics';
generalinfo.InstitutionAddress                      = 'Strasse des 17. Juni 135, 10623, Berlin, Germany';
generalinfo.TaskDescription                         = 'Participants walked repeatedly on a straight path.';
generalinfo.task                                    = bemobil_config.bids_tasklabel;  

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
subjectData.fields.age.Description      = 'age of the participant'; 
subjectData.fields.age.Unit             = 'years'; 
subjectData.fields.sex.Description      = 'sex of the participant'; 
subjectData.fields.sex.Levels.M         = 'male'; 
subjectData.fields.sex.Levels.F         = 'female'; 
subjectData.fields.group.Description    = 'experiment group';
subjectData.fields.group.Levels.young   = 'younger participants under 40';
subjectData.fields.group.Levels.old     = 'older participants over 65';
subjectData.fields.handedness.Description    = 'handedness of the participant';
subjectData.fields.handedness.Levels.R       = 'right-handed';
subjectData.fields.handedness.Levels.L       = 'left-handed';

% names of the columns - 'nr' column is just the numerical IDs of subjects
%                         do not change the name of this column 
subjectData.cols = {'nr',   'age',  'sex',  'group',    'handedness'};
subjectData.data = {64,     71,     'F',    'old',      'R' ; ...
                    66,     67,     'M',    'old',      'R' ; ...
                    76,     34,     'M',    'young',    'R' ; ...
                    78,     33,     'M',    'young',    'R' };


bemobil_xdf2bids(bemobil_config, generalinfo, subjectData)