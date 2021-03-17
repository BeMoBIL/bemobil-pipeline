% minimal version only for importing 

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
bemobil_config.bids_data_folder         = '1_BIDS-data_defaults\';
bemobil_config.bids_eegkeyword          = 'BrainVision';                  % marker streams also contain these strings. However, only the continuous stream is imported
bemobil_config.bids_tasklabel           = 'walking';

% custom function names - customization highly recommeded
bemobil_config.bids_motionconvert_custom    = 'bids_motionconvert_mobiworkshop';
bemobil_config.bids_parsemarkers_custom     = [];

% numerical IDs 
%--------------------------------------------------------------------------
numericalIDs                            = [64,66,76,78]; 

bemobil_xdf2bids(bemobil_config, numericalIDs) 