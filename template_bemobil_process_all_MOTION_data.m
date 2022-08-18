
% initialize EEGLAB 
if ~exist('ALLEEG','var')
	eeglab;
end

% load configuration 
template_bemobil_config;

% enter all subjects to process here (you can split it up in more MATLAB instances if you have more CPU power and RAM)
subjects = [1 2];

% set to 1 if all files should be computed, independently of whether they are present on disk or not
force_recompute = 0; 

%% processing loop

for subject = subjects
    
    %% prepare filepaths and check if already done
    
    disp(['Subject #' num2str(subject)]);
    
    STUDY = []; ALLEEG = [];  CURRENTSET=[]; MOTION_processed = [];
    
    input_filepath = [bemobil_config.study_folder bemobil_config.raw_EEGLAB_data_folder bemobil_config.filename_prefix num2str(subject)];
    output_filepath = [bemobil_config.study_folder bemobil_config.single_subject_motion_folder bemobil_config.filename_prefix num2str(subject)];
    
    try
        % load completely processed file
        MOTION_processed = pop_loadset('filename', [ bemobil_config.filename_prefix num2str(subject)...
            '_' bemobil_config.processed_motion_filename], 'filepath', output_filepath);
    catch
        disp('...failed. Computing now.')
    end
    
    if ~force_recompute && exist('MOTION_processed','var') && ~isempty(MOTION_processed)
        clear MOTION_processed
        disp('Subject is completely preprocessed already.')
        continue
    end
   	
	%% load data in EEGLAB .set structure
    % make sure the data is stored in double precision, large datafiles are supported, no memory mapped objects are
    % used but data is processed locally, and two files are used for storing sets (.set and .fdt)
	try 
        pop_editoptions('option_saveversion6', 0, 'option_single', 0, 'option_memmapdata', 0, 'option_savetwofiles', 1, 'option_storedisk', 0);
    catch
        warning('Could NOT edit EEGLAB memory options!!'); 
    end
    
    % load files that were created from xdf to BIDS to EEGLAB
    EEG_motion = pop_loadset('filename',[ bemobil_config.filename_prefix num2str(subject) '_' bemobil_config.raw_motion_filename],'filepath',input_filepath);
    
    %% individual EEG processing to remove non-exp segments
    % if you removed segments for the EEG processing, they have to be removed here, too, to keep the files synchronized
    
    % this example removes everything before the first and after the last event with a buffer of 1 second
    
    allevents = {EEG_motion.event.type}';
    
    removeindices = [0 EEG_motion.event(1).latency-EEG_motion.srate]; % remove from start to first event
    
    % add more removeIndices here for pauses or itnerruptions of the experiment if they have markers or you know their
    % indices in the data...
    
    removeindices(end+1,:) = [EEG_motion.event(end).latency+EEG_motion.srate EEG_motion.pnts]; % remove from last event to the end
    
    % filter for plot
    EEG_plot = pop_eegfiltnew(EEG_motion, 'locutoff',0.5,'plotfreqz',0);
    
    % plot
    fig1 = figure; set(gcf,'Color','w','InvertHardCopy','off', 'units','normalized','outerposition',[0 0 1 1])
    plot(normalize(EEG_plot.data') + [1:10:10*EEG_plot.nbchan], 'color', [78 165 216]/255)
    yticks([])
    
    xlim([0 EEG_motion.pnts])
    ylim([-10 10*EEG_plot.nbchan+10])
    
    hold on
    
    % plot lines for valid times
    
    for i = 1:size(removeindices,1)
        plot([removeindices(i,1) removeindices(i,1)],ylim,'r')
        plot([removeindices(i,2) removeindices(i,2)],ylim,'g')
    end
    
    % save plot
    print(gcf,fullfile(input_filepath,[bemobil_config.filename_prefix num2str(subject) '_raw-full_MOTION.png']),'-dpng')
    close

    % reject
    EEG_motion = eeg_eegrej(EEG_motion, removeindices);   
    
    %% processing wrapper for basic motion processing
  
    bemobil_process_all_motion(ALLEEG, EEG_motion, CURRENTSET, subject, bemobil_config, force_recompute);
    
end

subjects
subject

disp('PROCESSING DONE! YOU CAN CLOSE THE WINDOW NOW!')
