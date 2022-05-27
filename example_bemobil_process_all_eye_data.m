
% initialize EEGLAB 
if ~exist('ALLEEG','var')
	eeglab;
end

% load configuration 
example_bemobil_config;

% enter all subjects to process here (you can split it up in more MATLAB instances if you have more CPU power and RAM)
subjects = 1;

% set to 1 if all files should be computed, independently of whether they are present on disk or not
force_recompute = 0; 

%% additional parameters which are not based on the defaults but necessary in this case

bemobil_config.eye_raw_filename = 'merged_PHYSIO.set';
bemobil_config.eye_clean_filename = 'eye_clean.set';
bemobil_config.single_subject_eye_folder = ['7_single-subject-EYE-analysis' filesep];

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
    
    %% load data that is provided by the BIDS importer
    % make sure the data is stored in double precision, large datafiles are supported, and no memory mapped objects are
    % used but data is processed locally
	try 
        pop_editoptions( 'option_saveversion6', 0, 'option_single', 0, 'option_memmapdata', 0);
    catch
        warning('Could NOT edit EEGLAB memory options!!'); 
    end
    
    % load files that were created from xdf to BIDS to EEGLAB
    EYE = pop_loadset('filename',[ bemobil_config.filename_prefix num2str(subject) '_' bemobil_config.eye_raw_filename],'filepath',input_filepath);
    
    %% individual EEG processing to remove non-exp segments
    % if you removed segments for the EEG processing, they have to be removed here, too to keep the files synchronized
    
    %% processing eye 
    % the final EYE_clean data will contain additional blink events that can be easily transferred to other data files
    % if they are synchronized (which they should be when imported via the bemobil importing from xdf)
    
    idx_pupil = 7;
    idx_x = 1;
    idx_y = 2;
    
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
