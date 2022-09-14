% bemobil_process_all_EEG_preprocessing - wrapper function that incorporates all necessary processing steps from the basic EEG
% struct (e.g. all blocks merged together, nothing else done before except resampling) up to the preprocessed dataset
% which has line noise removed, channels interpolated, average reference, and relevant information stored in the EEG
% struct. Also plots several analytics plots along the way which are stored on disk alongside their respective files.
%
% Usage:
%   >>  [ALLEEG, EEG_interp_avRef, CURRENTSET] = bemobil_process_all_EEG_preprocessing(subject, bemobil_config, ALLEEG, EEG_merged,...
%     CURRENTSET, force_recompute)
%
% Inputs:
%   subject                   - subject number of the current subject (necessary for filepaths and storage)
%   bemobil_config            - configuration struct with all necessary information. See EEG_processing_example file
%                                that comes with this function!
%   ALLEEG                    - complete EEGLAB data set structure
%   EEG_to_process            - EEGLAB EEG structure that should be processed. Best to have all blocks merged into one
%                                file.
%   CURRENTSET                - index of current EEGLAB EEG structure within ALLEEG
%   force_recompute           - OPTIONAL force recomputation even if processed file is already present, default = 0
%
% Outputs:
%   ALLEEG                    - complete EEGLAB data set structure
%   EEG_interp_avRef          - processed EEGLAB EEG structure
%   Currentset                - index of current EEGLAB EEG structure within ALLEEG
%
%   .set data file of current EEGLAB EEG structure stored on disk (OPTIONALLY)
%
% See also:
%   EEGLAB, bemobil_process_EEG_basics, bemobil_detect_bad_channels, bemobil_interp_avref
%
% Authors: Marius Klug, 2021

function [ALLEEG, EEG_preprocessed, CURRENTSET] = bemobil_process_all_EEG_preprocessing(subject, bemobil_config, ALLEEG, EEG_to_process,...
    CURRENTSET, force_recompute)

% check config
bemobil_config = bemobil_check_config(bemobil_config);

%% basic setup

% make sure the data is stored in double precision, large datafiles are supported, no memory mapped objects are
% used but data is processed locally, and two files are used for storing sets (.set and .fdt)
try 
    pop_editoptions('option_saveversion6', 0, 'option_single', 0, 'option_memmapdata', 0, 'option_savetwofiles', 1, 'option_storedisk', 0);
catch
    warning('Could NOT edit EEGLAB memory options!!'); 
end

if ~exist('force_recompute','var') || isempty(force_recompute)
    force_recompute = 0;
end

disp(['Subject #' num2str(subject)]);


output_filepath = fullfile(bemobil_config.study_folder, bemobil_config.EEG_preprocessing_data_folder,...
    [bemobil_config.filename_prefix num2str(subject)]);
mkdir(output_filepath)


% check if the whole script has been running already
if ~force_recompute
    try
        EEG = pop_loadset('filename', [bemobil_config.filename_prefix num2str(subject) '_'...
            bemobil_config.preprocessed_filename], 'filepath', output_filepath);
        [ALLEEG, EEG_preprocessed, CURRENTSET] = pop_newset(ALLEEG, EEG, 0,'study',0);
        
        warning('Old interpolated file already existed, using that file!')
        
        return
    catch
        disp('...failed. Computing now.')
    end
end

%% EEG basics - resample, chanlocs, zapline, add ref with zeros, declare channel type
if ~force_recompute
    try
        
        EEG_basic = pop_loadset('filename', [bemobil_config.filename_prefix num2str(subject) '_'...
            bemobil_config.basic_prepared_filename], 'filepath', output_filepath);
        
        warning('Old basic prepared file already existed, using that file!')
        
    catch
        disp('...failed. Computing now.')
    end
end

if ~exist('EEG_basic','var')
    
    %% plot raw data at 6 different times throughout the dataset
    % 0.5hz filter just for plotting, so remove DC offset
    
    disp('Filtering data only for plotting!')
    EEG = pop_eegfiltnew(EEG_to_process, 'locutoff',0.5);
    
    %% create plots
    
    plotfigure = figure('color','w');
    set(plotfigure, 'Position', get(0,'screensize'))
    ax1 = subplot(231);
    ax2 = subplot(232);
    ax3 = subplot(233);
    ax4 = subplot(234);
    ax5 = subplot(235);
    ax6 = subplot(236);
    
    starttime = EEG.times(end)/7*1;
    vis_artifacts(EEG,EEG,'show_events',1,'time_subset',...
        round([starttime starttime+10000]/1000)); % plot 10s at the first quarter
    drawnow
    axeshandle = gca;
    fighandle = gcf;
    axcp = copyobj(axeshandle, plotfigure);
    set(axcp,'Position',get(ax1,'position'));
    axcp.XTickLabel = [0:10]+round(starttime/1000);
    axcp.YTick=[];
    axcp.Title.String = ['Raw data section 1 of ' num2str(round(EEG.times(end)/1000)) 's'];
    axcp.XLabel.String = 'seconds';
    drawnow
    delete(ax1);
    close(fighandle)
    
    starttime = EEG.times(end)/7*2;
    vis_artifacts(EEG,EEG,'show_events',1,'time_subset',...
        round([starttime starttime+10000]/1000)); % plot 10s at the first quarter
    drawnow
    axeshandle = gca;
    fighandle = gcf;
    axcp = copyobj(axeshandle, plotfigure);
    set(axcp,'Position',get(ax2,'position'));
    axcp.XTickLabel = [0:10]+round(starttime/1000);
    axcp.YTick=[];
    axcp.Title.String = ['Raw data section 2 of ' num2str(round(EEG.times(end)/1000)) 's'];
    axcp.XLabel.String = 'seconds';
    drawnow
    delete(ax2);
    close(fighandle)
    
    starttime = EEG.times(end)/7*3;
    vis_artifacts(EEG,EEG,'show_events',1,'time_subset',...
        round([starttime starttime+10000]/1000)); % plot 10s at the first quarter
    drawnow
    axeshandle = gca;
    fighandle = gcf;
    axcp = copyobj(axeshandle, plotfigure);
    set(axcp,'Position',get(ax3,'position'));
    axcp.XTickLabel = [0:10]+round(starttime/1000);
    axcp.YTick=[];
    axcp.Title.String = ['Raw data section 3 of ' num2str(round(EEG.times(end)/1000)) 's'];
    axcp.XLabel.String = 'seconds';
    drawnow
    delete(ax3);
    close(fighandle)
    
    starttime = EEG.times(end)/7*4;
    vis_artifacts(EEG,EEG,'show_events',1,'time_subset',...
        round([starttime starttime+10000]/1000)); % plot 10s at the first quarter
    drawnow
    axeshandle = gca;
    fighandle = gcf;
    axcp = copyobj(axeshandle, plotfigure);
    set(axcp,'Position',get(ax4,'position'));
    axcp.XTickLabel = [0:10]+round(starttime/1000);
    axcp.YTick=[];
    axcp.Title.String = ['Raw data section 4 of ' num2str(round(EEG.times(end)/1000)) 's'];
    axcp.XLabel.String = 'seconds';
    drawnow
    delete(ax4);
    close(fighandle)
    
    starttime = EEG.times(end)/7*5;
    vis_artifacts(EEG,EEG,'show_events',1,'time_subset',...
        round([starttime starttime+10000]/1000)); % plot 10s at the first quarter
    drawnow
    axeshandle = gca;
    fighandle = gcf;
    axcp = copyobj(axeshandle, plotfigure);
    set(axcp,'Position',get(ax5,'position'));
    axcp.XTickLabel = [0:10]+round(starttime/1000);
    axcp.YTick=[];
    axcp.Title.String = ['Raw data section 5 of ' num2str(round(EEG.times(end)/1000)) 's'];
    axcp.XLabel.String = 'seconds';
    drawnow
    delete(ax5);
    close(fighandle)
    
    starttime = EEG.times(end)/7*6;
    vis_artifacts(EEG,EEG,'show_events',1,'time_subset',...
        round([starttime starttime+10000]/1000)); % plot 10s at the first quarter
    drawnow
    axeshandle = gca;
    fighandle = gcf;
    axcp = copyobj(axeshandle, plotfigure);
    set(axcp,'Position',get(ax6,'position'));
    axcp.XTickLabel = [0:10]+round(starttime/1000);
    axcp.YTick=[];
    axcp.Title.String = ['Raw data section 6 of ' num2str(round(EEG.times(end)/1000)) 's'];
    axcp.XLabel.String = 'seconds';
    drawnow
    delete(ax6);
    close(fighandle)
    
    
    %% save plot
    
    savefig(plotfigure,fullfile(output_filepath,[bemobil_config.filename_prefix num2str(subject) '_raw.fig']))
    print(plotfigure,fullfile(output_filepath,[bemobil_config.filename_prefix num2str(subject) '_raw.png']),'-dpng')
    close
    
    if ~isempty(bemobil_config.channel_locations_filename)
        channel_locations_filepath = fullfile(bemobil_config.study_folder, bemobil_config.source_data_folder,...
            [bemobil_config.filename_prefix num2str(subject)], [bemobil_config.filename_prefix num2str(subject) '_'...
            bemobil_config.channel_locations_filename]);
    else
        channel_locations_filepath = [];
    end
    
    % preprocessing: enter chanlocs, remove unused channels, declare EOG, resample
    [ALLEEG, EEG_basic, CURRENTSET] = bemobil_process_EEG_basics(ALLEEG, EEG_to_process, CURRENTSET, channel_locations_filepath,...
        bemobil_config.channels_to_remove, bemobil_config.eog_channels, bemobil_config.resample_freq,...
        [bemobil_config.filename_prefix num2str(subject) '_' bemobil_config.basic_prepared_filename], output_filepath,...
        bemobil_config.rename_channels, bemobil_config.ref_channel, bemobil_config.zaplineConfig);
    
    disp('Preprocessing done!')
    
end

% save RAM
clear EEG_to_process

%% detect bad channels
[chans_to_interp, chan_detected_fraction_threshold, detected_bad_channels, rejected_chan_plot_handle, detection_plot_handle] = bemobil_detect_bad_channels(EEG_basic, ALLEEG, CURRENTSET,...
    bemobil_config.chancorr_crit,bemobil_config.chan_max_broken_time, bemobil_config.chan_detect_num_iter,...
    bemobil_config.chan_detected_fraction_threshold,bemobil_config.num_chan_rej_max_target, bemobil_config.flatline_crit,bemobil_config.line_noise_crit);

if length(chans_to_interp) > EEG_basic.nbchan/5
    warndlg(['In subject ' num2str(subject) ', ' num2str(length(chans_to_interp)) ' of ' num2str(EEG_basic.nbchan)...
        ' channels were rejected, which is more than 1/5th!'])
end

EEG_basic.etc.channel_rejection.detection_threshold = chan_detected_fraction_threshold;
EEG_basic.etc.channel_rejection.bad_channel_detection = detected_bad_channels;

%% save fig of bad channels

savefig(rejected_chan_plot_handle,fullfile(output_filepath,[bemobil_config.filename_prefix num2str(subject) '_bad_channels.fig']))
print(rejected_chan_plot_handle,fullfile(output_filepath,[bemobil_config.filename_prefix num2str(subject) '_bad_channels.png']),'-dpng')
close(rejected_chan_plot_handle)

savefig(detection_plot_handle,fullfile(output_filepath,[bemobil_config.filename_prefix num2str(subject) '_bad_channels_detection.fig']))
print(detection_plot_handle,fullfile(output_filepath,[bemobil_config.filename_prefix num2str(subject) '_bad_channels_detection.png']),'-dpng')
close(detection_plot_handle)

%% do the actual interpolation and average referencing (reference is not considering EOGs)

disp('Interpolating bad channels and compute final average reference, ignoring EOG channels...')
[ALLEEG, EEG_preprocessed, CURRENTSET] = bemobil_interp_avref( EEG_basic , ALLEEG, CURRENTSET, chans_to_interp,...
    [bemobil_config.filename_prefix num2str(subject) '_' bemobil_config.preprocessed_filename], output_filepath);

%% plot interpolated filtered, for analytics

disp('Filtering data only for plotting!')
EEG = pop_eegfiltnew(EEG_preprocessed, 'locutoff',0.5);

%%

plotfigure = figure('color','w');
set(plotfigure, 'Position', get(0,'screensize'))
ax1 = subplot(231);
ax2 = subplot(232);
ax3 = subplot(233);
ax4 = subplot(234);
ax5 = subplot(235);
ax6 = subplot(236);

starttime = EEG.times(end)/7*1;
vis_artifacts(EEG,EEG,'show_events',1,'time_subset',...
    round([starttime starttime+10000]/1000)); % plot 10s at the first quarter
drawnow
axeshandle = gca;
fighandle = gcf;
axcp = copyobj(axeshandle, plotfigure);
set(axcp,'Position',get(ax1,'position'));
axcp.XTickLabel = [0:10]+round(starttime/1000);
axcp.YTick=[];
axcp.Title.String = ['Interpolated channels data section 1 of ' num2str(round(EEG.times(end)/1000)) 's'];
axcp.XLabel.String = 'seconds';
drawnow
delete(ax1);
close(fighandle)

starttime = EEG.times(end)/7*2;
vis_artifacts(EEG,EEG,'show_events',1,'time_subset',...
    round([starttime starttime+10000]/1000)); % plot 10s at the first quarter
drawnow
axeshandle = gca;
fighandle = gcf;
axcp = copyobj(axeshandle, plotfigure);
set(axcp,'Position',get(ax2,'position'));
axcp.XTickLabel = [0:10]+round(starttime/1000);
axcp.YTick=[];
axcp.Title.String = ['Interpolated channels data section 2 of ' num2str(round(EEG.times(end)/1000)) 's'];
axcp.XLabel.String = 'seconds';
drawnow
delete(ax2);
close(fighandle)

starttime = EEG.times(end)/7*3;
vis_artifacts(EEG,EEG,'show_events',1,'time_subset',...
    round([starttime starttime+10000]/1000)); % plot 10s at the first quarter
drawnow
axeshandle = gca;
fighandle = gcf;
axcp = copyobj(axeshandle, plotfigure);
set(axcp,'Position',get(ax3,'position'));
axcp.XTickLabel = [0:10]+round(starttime/1000);
axcp.YTick=[];
axcp.Title.String = ['Interpolated channels data section 3 of ' num2str(round(EEG.times(end)/1000)) 's'];
axcp.XLabel.String = 'seconds';
drawnow
delete(ax3);
close(fighandle)

starttime = EEG.times(end)/7*4;
vis_artifacts(EEG,EEG,'show_events',1,'time_subset',...
    round([starttime starttime+10000]/1000)); % plot 10s at the first quarter
drawnow
axeshandle = gca;
fighandle = gcf;
axcp = copyobj(axeshandle, plotfigure);
set(axcp,'Position',get(ax4,'position'));
axcp.XTickLabel = [0:10]+round(starttime/1000);
axcp.YTick=[];
axcp.Title.String = ['Interpolated channels data section 4 of ' num2str(round(EEG.times(end)/1000)) 's'];
axcp.XLabel.String = 'seconds';
drawnow
delete(ax4);
close(fighandle)

starttime = EEG.times(end)/7*5;
vis_artifacts(EEG,EEG,'show_events',1,'time_subset',...
    round([starttime starttime+10000]/1000)); % plot 10s at the first quarter
drawnow
axeshandle = gca;
fighandle = gcf;
axcp = copyobj(axeshandle, plotfigure);
set(axcp,'Position',get(ax5,'position'));
axcp.XTickLabel = [0:10]+round(starttime/1000);
axcp.YTick=[];
axcp.Title.String = ['Interpolated channels data section 5 of ' num2str(round(EEG.times(end)/1000)) 's'];
axcp.XLabel.String = 'seconds';
drawnow
delete(ax5);
close(fighandle)

starttime = EEG.times(end)/7*6;
vis_artifacts(EEG,EEG,'show_events',1,'time_subset',...
    round([starttime starttime+10000]/1000)); % plot 10s at the first quarter
drawnow
axeshandle = gca;
fighandle = gcf;
axcp = copyobj(axeshandle, plotfigure);
set(axcp,'Position',get(ax6,'position'));
axcp.XTickLabel = [0:10]+round(starttime/1000);
axcp.YTick=[];
axcp.Title.String = ['Interpolated channels data section 6 of ' num2str(round(EEG.times(end)/1000)) 's'];
axcp.XLabel.String = 'seconds';
drawnow
delete(ax6);
close(fighandle)


%% save plot

savefig(plotfigure,fullfile(output_filepath,[bemobil_config.filename_prefix num2str(subject) '_interpolated_channels.fig']))
print(plotfigure,fullfile(output_filepath,[bemobil_config.filename_prefix num2str(subject) '_interpolated_channels.png']),'-dpng')
close

disp('All basic EEG processing done.')