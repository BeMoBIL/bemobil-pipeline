
% initialize EEGLAB 
if ~exist('ALLEEG','var')
	eeglab;
end

% load configuration 
example_bemobil_config;

% enter all subjects to process here (you can split it up in more MATLAB instances if you have more CPU power and RAM)
subjects = [64 66 76 78];

% set to 1 if all files should be computed, independently of whether they are present on disk or not
force_recompute = 0; 

%% processing loop

for subject = subjects
    
    %% prepare filepaths and check if already done
    
    disp(['Subject #' num2str(subject)]);
    
    STUDY = []; ALLEEG = [];  CURRENTSET=[]; motion_processed = [];
    
    input_filepath = [bemobil_config.study_folder bemobil_config.raw_EEGLAB_data_folder bemobil_config.filename_prefix num2str(subject)];
    output_filepath = [bemobil_config.study_folder bemobil_config.motion_analysis_folder bemobil_config.filename_prefix num2str(subject)];
    
    try
        % load completely processed file
        motion_processed = pop_loadset('filename', [ bemobil_config.filename_prefix num2str(subject)...
            '_' bemobil_config.processed_motion_filename], 'filepath', output_filepath);
    catch
        disp('...failed.')
    end
    
    if ~force_recompute && exist('motion_processed','var') && ~isempty(motion_processed)
        clear motion_processed
        disp('Subject is completely preprocessed already.')
        continue
    end
    
    %% load data that is provided by the BIDS importer
    
    % make sure the data is stored in double precision, large datafiles are supported, no memory mapped objects are
    % used but data is processed locally, and two files are used for storing sets (.set and .fdt)
    try
        pop_editoptions('option_saveversion6', 0, 'option_single', 0, 'option_memmapdata', 0, 'option_savetwofiles', 1, 'option_storedisk', 0);
    catch
        warning('Could NOT edit EEGLAB memory options!!');
    end
    
    % load files that were created from xdf to BIDS to EEGLAB
    EEG = pop_loadset('filename',[ bemobil_config.filename_prefix num2str(subject) '_RAW_MOTION_DATA_FILENAME.set'],'filepath',input_filepath);
    
    %% individual EEG processing to remove non-exp segments
    % if you removed segments for the EEG processing, they have to be removed here, too to keep the files synchronized
    
    %% processing wrapper for basic motion processing
  
    bemobil_process_all_motion(ALLEEG, EEG, CURRENTSET, subject, bemobil_config, force_recompute);
    
end

subjects
subject

disp('PROCESSING DONE! YOU CAN CLOSE THE WINDOW NOW!')
