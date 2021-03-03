% example script for modality agnostic metadata
% metadata are dataset-specific
%--------------------------------------------------------------------------
% initialize general config struct shared across modalities
cfg = [];

% root directory (where you want your bids data to be saved)
cfg.bidsroot                                = 'E:\Project_BIDS\example_dataset_MWM\rawdata_bids'; 

% required for dataset_description.json
cfg.dataset_description.Name                = 'Data set collected at BeMoBIL';
cfg.dataset_description.BIDSVersion         = 'unofficial extension';

% optional for dataset_description.json
cfg.dataset_description.License             = 'n/a';
cfg.dataset_description.Authors             = '';
cfg.dataset_description.Acknowledgements    = 'Acknowledgements here';
cfg.dataset_description.Funding             = 'n/a';
cfg.dataset_description.ReferencesAndLinks  = 'n/a';
cfg.dataset_description.DatasetDOI          = 'n/a';

% general information shared across modality specific json files 
cfg.InstitutionName                         = 'Technische Universitaet zu Berlin';
cfg.InstitutionalDepartmentName             = 'Biological Psychology and Neuroergonomics';
cfg.InstitutionAddress                      = 'Strasse des 17. Juni 135, 10623, Berlin, Germany';
cfg.TaskDescription                         = 'Here description of your task.';
cfg.task                                    = 'taskname';                   % this is the string to describe your task, to go in file names