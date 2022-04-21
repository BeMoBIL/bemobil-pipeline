% initialize EEGLAB 
if ~exist('ALLCOM','var')
	eeglab;
end

% load configuration 
example_bemobil_config;

% enter all subjects to process here (you can split it up in more MATLAB instances if you have more CPU power and RAM)
subjects = 1:2; 

% set to 1 if all files should be computed, independently of whether they are present on disk or not
force_recompute = 0; 

%% Import 

% not in this example 

%% processing loop

for subject = subjects
    
    %% prepare filepaths and check if already done
    
	disp(['Subject #' num2str(subject)]);
    
	STUDY = []; CURRENTSTUDY = 0; ALLEEG = [];  CURRENTSET=[]; EEG=[]; EEG_interp_avref = []; EEG_single_subject_final = [];
	
	input_filepath = [bemobil_config.study_folder bemobil_config.raw_EEGLAB_data_folder bemobil_config.filename_prefix num2str(subject)];
	output_filepath = [bemobil_config.study_folder bemobil_config.single_subject_analysis_folder bemobil_config.filename_prefix num2str(subject)];
	
	try
		% load completely processed file
		EEG_single_subject_final = pop_loadset('filename', [ bemobil_config.filename_prefix num2str(subject)...
			'_' bemobil_config.single_subject_cleaned_ICA_filename], 'filepath', output_filepath);
    catch
        disp('...failed. Computing now.')
	end
	
	if ~force_recompute && exist('EEG_single_subject_final','var') && ~isempty(EEG_single_subject_final)
		clear EEG_single_subject_final
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
    EEG = pop_loadset('filename',[ bemobil_config.filename_prefix num2str(subject) '_' bemobil_config.merged_filename],'filepath',input_filepath);
    
    %% individual EEG processing to remove non-exp segments
    % it is stongly recommended to remove these because they may contain strong artifacts that confuse channel detection
    % and AMICA
    
    % this example removes everything before the first and after the last event, plus is searches for the pause keyword
    % to remove breaks with a buffer of 1 second
    
    allevents = {EEG.event.type}';
    
    idx_startpause = find(contains(lower(allevents),'start:pause'));
    idx_stoppause = find(contains(lower(allevents),'stop:pause'));
    
    if (length(idx_startpause) ~= length(idx_stoppause))
        if (length(idx_startpause) == length(idx_stoppause)+1) && (idx_startpause(end) > idx_stoppause(end))
            warning('Last stop pause missing, removing last start pause!')
            idx_startpause(end) = [];
        end
    end
    
    removeindices = [0 EEG.event(1).latency-EEG.srate];
    
    for i = 1:length(idx_startpause)
        removeindices(end+1,1) = round(EEG.event(idx_startpause(i)).latency)+EEG.srate;
        removeindices(end,2) = round(EEG.event(idx_stoppause(i)).latency)-EEG.srate;
    end
    
    removeindices(end+1,:) = [EEG.event(end).latency+EEG.srate EEG.pnts];
    
    % filter for plot
    EEG_plot = pop_eegfiltnew(EEG, 'locutoff',0.5,'plotfreqz',0);
    
    % plot
    fig1 = figure; set(gcf,'Color','w','InvertHardCopy','off', 'units','normalized','outerposition',[0 0 1 1])

    % basic chan reject for plot
    chanmaxes = max(EEG_plot.data,[],2);
    EEG_plot = pop_select( EEG_plot, 'nochannel',find(chanmaxes>median(chanmaxes)+1.4826 * 3* mad(chanmaxes)));
    chanmins = min(EEG_plot.data,[],2);
    EEG_plot = pop_select( EEG_plot, 'nochannel',find(chanmins<median(chanmins)-1.4826 * 3* mad(chanmins)));
    
    plot(EEG_plot.data' + linspace(0,20000,EEG_plot.nbchan), 'color', [78 165 216]/255)
    xlim([0 EEG.pnts])
    ylim([-1000 21000])
    
    hold on
    
    % plot lines for valid times
    
    for i = 1:size(removeindices,1)
        plot([removeindices(i,1) removeindices(i,1)],[-1000 21000],'r')
        plot([removeindices(i,2) removeindices(i,2)],[-1000 21000],'g')
    end
    
    % save plot
    print(gcf,fullfile(input_filepath,[bemobil_config.filename_prefix num2str(subject) '_raw-full.png']),'-dpng')
    close

    % reject
    EEG = eeg_eegrej(EEG, removeindices);   
    
    %% processing wrappers for basic processing and AMICA
    
    % do basic preprocessing, line noise removal, and channel interpolation
	[ALLEEG, EEG_preprocessed, CURRENTSET] = bemobil_process_all_EEG_preprocessing(subject, bemobil_config, ALLEEG, EEG, CURRENTSET, force_recompute);

    % start the processing pipeline for AMICA
	bemobil_process_all_AMICA(ALLEEG, EEG_preprocessed, CURRENTSET, subject, bemobil_config, force_recompute);

end

subjects
subject

disp('PROCESSING DONE! YOU CAN CLOSE THE WINDOW NOW!')
