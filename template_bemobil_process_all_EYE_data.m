
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
    
	STUDY = []; CURRENTSTUDY = 0; ALLEEG = [];  CURRENTSET=[]; EYE=[]; EYE_clean = [];
	
	input_filepath = [bemobil_config.study_folder bemobil_config.raw_EEGLAB_data_folder bemobil_config.filename_prefix num2str(subject)];
	output_filepath = [bemobil_config.study_folder bemobil_config.single_subject_eye_folder bemobil_config.filename_prefix num2str(subject)];
	
	try
		% load completely processed file
		EYE_clean = pop_loadset('filename', [ bemobil_config.filename_prefix num2str(subject)...
			'_' bemobil_config.eye_clean_filename], 'filepath', output_filepath);
    catch
        disp('...failed. Computing now.')
	end
	
	if ~force_recompute && exist('EYE_clean','var') && ~isempty(EYE_clean)
		clear EYE_clean
		disp('Subject is completely preprocessed already.')
		continue
    end
    
    %% load data in EEGLAB .set structure
    % make sure the data is stored in double precision, large datafiles are supported, no memory mapped objects are
    % used but data is processed locally, and two files are used for storing sets (.set and .fdt)
	try 
        pop_editoptions( 'option_saveversion6', 0, 'option_single', 0, 'option_memmapdata', 0);
    catch
        warning('Could NOT edit EEGLAB memory options!!'); 
    end
    
    % load files that were created from xdf to BIDS to EEGLAB
    EYE = pop_loadset('filename',[ bemobil_config.filename_prefix num2str(subject) '_' bemobil_config.eye_raw_filename],'filepath',input_filepath);
    
    %% individual EEG processing to remove non-exp segments
    % if you removed segments for the EEG processing, they have to be removed here, too to keep the files synchronized
    
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
    print(gcf,fullfile(input_filepath,[bemobil_config.filename_prefix num2str(subject) '_raw-full_EYE.png']),'-dpng')
    close

    % reject
    EEG_motion = eeg_eegrej(EEG_motion, removeindices);   
    
    %% processing eye 
    % the final EYE_clean data will contain additional blink events that can be easily transferred to other data files
    % if they are synchronized (which they should be when imported via the bemobil importing from xdf)
    
    idx_pupil = [];
    idx_x = [];
    idx_y = [];
    
    [EYE_clean, plothandles] = bemobil_clean_eye(EYE,idx_pupil);
    
    % add figure for x-y plot
    xyplothandle = figure('color','w','position',[50 100 1500 800]);
    subplot(121)
    plot(EYE.data(idx_x,:),EYE.data(idx_y,:))
    title('Raw')
    xlabel('Screen X both')
    ylabel('Screen Y both')
    subplot(122)
    plot(EYE_clean.data(idx_x,:),EYE_clean.data(idx_y,:))
    title('Cleaned')
    xlabel('Screen X both')
    ylabel('Screen Y both')
    
    % save data and plots
    mkdir(output_filepath)
    pop_saveset(EYE_clean, 'filename', [ bemobil_config.filename_prefix num2str(subject)...
		'_' bemobil_config.eye_clean_filename], 'filepath', output_filepath);
    disp('done.')
    
    print(xyplothandle,fullfile(output_filepath,[bemobil_config.filename_prefix num2str(subject) '_xyplot.png']),'-dpng')
    savefig(xyplothandle,fullfile(output_filepath,[bemobil_config.filename_prefix num2str(subject) '_xyplot.fig']))
    close(xyplothandle)
    
    print(plothandles(1),fullfile(output_filepath,[bemobil_config.filename_prefix num2str(subject) '_raw.png']),'-dpng')
 
    print(plothandles(2),fullfile(output_filepath,[bemobil_config.filename_prefix num2str(subject) '_pupilradius.png']),'-dpng')
    savefig(plothandles(2),fullfile(output_filepath,[bemobil_config.filename_prefix num2str(subject) '_pupilradius.fig']))
    
    print(plothandles(3),fullfile(output_filepath,[bemobil_config.filename_prefix num2str(subject) '_clean.png']),'-dpng')
    
    for i=1:3
        close(plothandles(i))
    end
        
    disp('All figures and data saved!')
    
end

subjects
subject

disp('PROCESSING DONE! YOU CAN CLOSE THE WINDOW NOW!')
