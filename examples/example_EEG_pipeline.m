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
    
	STUDY = []; CURRENTSTUDY = 0; ALLEEG = [];  CURRENTSET=[]; EEG=[]; EEG_interp_avref = []; EEG_single_subject_final = [];
	
	input_filepath = [bemobil_config.study_folder bemobil_config.raw_EEGLAB_data_folder bemobil_config.filename_prefix num2str(subject)];
	output_filepath = [bemobil_config.study_folder bemobil_config.single_subject_analysis_folder bemobil_config.filename_prefix num2str(subject)];
	
	try
		% load completely processed file
		EEG_single_subject_final = pop_loadset('filename', [ bemobil_config.filename_prefix num2str(subject)...
			'_' erase(bemobil_config.preprocessed_and_ICA_filename,'.set') '_filtered.set'], 'filepath', output_filepath);
    catch
        disp('...failed. Computing now.')
	end
	
	if ~force_recompute && exist('EEG_single_subject_final','var') && ~isempty(EEG_single_subject_final)
		clear EEG_single_subject_final
		disp('Subject is completely preprocessed already.')
		continue
    end
	
	%% load data that is provided by the BIDS importer
    % make sure the data is stored in double precision, large datafiles are supported, and no memory mapped objects are
    % used but data is processed locally
	try 
        pop_editoptions( 'option_saveversion6', 0, 'option_single', 0, 'option_memmapdata', 0);
    catch
        warning('Could NOT edit EEGLAB memory options!!'); 
    end
    
    % load files that were created from xdf to BIDS to EEGLAB
    EEG = pop_loadset('filename',[ bemobil_config.filename_prefix num2str(subject) '_' bemobil_config.merged_filename],'filepath',input_filepath);
    
    %% individual EEG processing to remove non-exp segments
    
    allevents = {EEG.event.type}';

    startevents = {EEG.event(find(~cellfun(@isempty,strfind(allevents,'START')) & cellfun(@isempty,strfind(allevents,'TES')))).type}';
    startlatencies = [EEG.event(find(~cellfun(@isempty,strfind(allevents,'START')) & cellfun(@isempty,strfind(allevents,'TES')))).latency]';
    endevents = {EEG.event(find(~cellfun(@isempty,strfind(allevents,'END')) & cellfun(@isempty,strfind(allevents,'test')))).type}';
    endlatencies = [EEG.event(find(~cellfun(@isempty,strfind(allevents,'END')) & cellfun(@isempty,strfind(allevents,'test')))).latency]';

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

    latencies = [[1; endlatencies+EEG.srate] [startlatencies-EEG.srate; EEG.pnts]]; % remove segments but leave a buffer of 1 sec before and after the events for timefreq analysis

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
    
    for i = 1:size(latencies,1)
        plot([latencies(i,1) latencies(i,1)],[-1000 21000],'r')
        plot([latencies(i,2) latencies(i,2)],[-1000 21000],'g')
    end
    
    % save plot
    print(gcf,fullfile(input_filepath,[bemobil_config.filename_prefix num2str(subject) '_' erase(bemobil_config.merged_filename,'.set') '_raw-full.png']),'-dpng')
    close
    
    EEG = eeg_eegrej( EEG, latencies);
    
    %% processing wrappers for basic stuff and AMICA
    
    % do basic preprocessing, line noise removal, and channel interpolation
	[ALLEEG, EEG_preprocessed, CURRENTSET] = bemobil_process_all_EEG_preprocessing(subject, bemobil_config, ALLEEG, EEG, CURRENTSET, force_recompute);

    % start the processing pipeline for AMICA
	bemobil_process_all_AMICA(ALLEEG, EEG_preprocessed, CURRENTSET, subject, bemobil_config, force_recompute);

end

subjects
subject

bemobil_copy_plots_in_one(bemobil_config)

disp('PROCESSING DONE! YOU CAN CLOSE THE WINDOW NOW!')
