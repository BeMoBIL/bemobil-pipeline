% example script for BeMoBIL BIDS tools bids2set

% configuration 
bemobil_config.study_folder             = 'E:\Project_BIDS\example_dataset_MWM\';
bemobil_config.filename_prefix          = 'sub';
bemobil_config.bids_data_folder         = '1_BIDS-data\';
bemobil_config.raw_EEGLAB_data_folder   = '2_basic-EEGLAB\';
bemobil_config.filenames                = {'desktop' 'VR'}; 

% these streams should be processed as rigid body streams containing 3 dof position and 3 dof orientation data (e.g. derivatives and filters applied)
bemobil_config.bids_tasklabel           = 'VNE1';

bemobil_bids2set(bemobil_config)