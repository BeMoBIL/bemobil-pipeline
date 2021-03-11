% example script for BeMoBIL BIDS tools xdf2bids

% configuration 
numericalIDs                            = [81001, 81002, 82001, 82002];
bemobil_config.study_folder             = 'E:\Project_BIDS\example_dataset_MWM\';
bemobil_config.filename_prefix          = 'sub';
bemobil_config.raw_data_folder          = 'sourcedata\';
bemobil_config.bids_data_folder         = '1_BIDS-data\';
bemobil_config.raw_EEGLAB_data_folder   = '2_basic-EEGLAB\';
bemobil_config.filenames                = {'desktop' 'VR'}; 

% these streams should be processed as rigid body streams containing 3 dof position and 3 dof orientation data (e.g. derivatives and filters applied)
bemobil_config.rigidbody_streams        = {'PlayerTransform','RightFoot', 'LeftFoot', 'Torso'};
bemobil_config.bids_rbsessions          = logical([1,1,0,0,0;1,1,1,1,1]); 
bemobil_config.bids_eegkeyword          = {'BrainVision'};
bemobil_config.channel_locations_filename = 'VN_E1_eloc.elc'; 
bemobil_config.bids_tasklabel           = 'VNE1';

bemobil_xdf2bids(bemobil_config, numericalIDs)