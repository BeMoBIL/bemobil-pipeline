bidsFolder      = '...\bids-data'; % root directory (where you want your bids data to be saved)
setFolder       = '...\2_EEGlab-basic'; % where you want to save the resulting set files
sessionNames    = {'ses1', 'ses2'}; % names of your sessions. If only one session is included, write task name  

% enter all subjects to process here 
subjects = [1 2 3 4];


%%  loop over subjects 
for subject = subjects
    
    config                          = []; % reset for each loop
    config.bids_target_folder       = bidsFolder; % required
    config.subject                  = subject; % required
    config.set_folder             	= setFolder;
    config.session_names            = sessionNames; 
    config.other_data_types         = {'motion'}; % other data type than eeg                      
    bemobil_bids2set(config);
    
end

