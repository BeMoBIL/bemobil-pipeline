% close all; clear

example_bemobil_config_script;

% enter all subjects to process here (you can split it up in more MATLAB instances if you have more CPU power and RAM)
subjects = [64 66 76 78];

force_recompute = 0;

%% processing loop

if ~exist('ALLCOM','var')
    eeglab;
end

%%

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
    % make sure the data is stored in double precision, large datafiles are supported, and no memory mapped objects are
    % used but data is processed locally
    try
        pop_editoptions('option_saveversion6', 0, 'option_single', 0, 'option_memmapdata', 0, 'option_savetwofiles', 1, 'option_storedisk', 0);
    catch
        warning('Could NOT edit EEGLAB memory options!!');
    end
    
    % load files that were created from xdf to BIDS to EEGLAB
    MOTION = pop_loadset('filename',[ bemobil_config.filename_prefix num2str(subject) '_' bemobil_config.merged_motion_filename],'filepath',input_filepath);
    
    %% remove non-exp segments as in the EEG processing
    
    allevents = {MOTION.event.type}';
    
    startevents = {MOTION.event(find(~cellfun(@isempty,strfind(allevents,'START')) & cellfun(@isempty,strfind(allevents,'TES')))).type}';
    startlatencies = [MOTION.event(find(~cellfun(@isempty,strfind(allevents,'START')) & cellfun(@isempty,strfind(allevents,'TES')))).latency]';
    endevents = {MOTION.event(find(~cellfun(@isempty,strfind(allevents,'END')) & cellfun(@isempty,strfind(allevents,'test')))).type}';
    endlatencies = [MOTION.event(find(~cellfun(@isempty,strfind(allevents,'END')) & cellfun(@isempty,strfind(allevents,'test')))).latency]';
    
    switch subject
        case 66
            startlatencies = startlatencies([1:7 9:end]);
            startevents = startevents([1:7 9:end]);
        case 76 % this subject had some issues in the first baseline
            startevents = startevents([2: 9 12:end]);
            startlatencies = startlatencies([2: 9 12:end]);
            endevents = endevents([2: 9 12:end]);
            endlatencies = endlatencies([2: 9 12:end]);
    end
    
    t = table(startevents,...
        startlatencies,...
        endlatencies,...
        endevents,...
        endlatencies - startlatencies);
    
    latencies = [[1; endlatencies+MOTION.srate] [startlatencies-MOTION.srate; MOTION.pnts]]; % remove segments but leave a buffer of 1 sec before and after the events for timefreq analysis
    
    MOTION = eeg_eegrej( MOTION, latencies);
    
    %% processing wrapper for motion
    
    bemobil_process_all_motion(ALLEEG, MOTION, CURRENTSET, subject, bemobil_config, force_recompute);
    
end

subjects
subject

disp('PROCESSING DONE! YOU CAN CLOSE THE WINDOW NOW!')
