% bemobil_process_all_AMICA - wrapper function that incorporates all necessary processing steps from the preprocessed
% EEG struct (line noise removed, bad channels interpolated, average referenced) to AMICA computation, dipole fitting,
% and artifact IC cleaning. AMICA autorejection is supported and recommended, all information including the rejected
% time points is stored in the final EEG set. A processing config struct is necessary. For an example please see the
% EEG_processing_example script! Plots several analytics figures that are stored alongside their respective datasets.
%
% The intermediate files are stored on the disk.
%
% Usage:
%   >>  [ALLEEG, EEG_preprocessed_and_ICA, CURRENTSET] = bemobil_process_all_AMICA(ALLEEG, EEG_interp_avRef, CURRENTSET,...
%     subject, bemobil_config, force_recompute)
% 
% Inputs:
%   ALLEEG                    - complete EEGLAB data set structure
%   EEG_to_process            - EEGLAB EEG structure that should be processed. Best to have all blocks merged into one
%                                file.
%   CURRENTSET                - index of current EEGLAB EEG structure within ALLEEG
%   subject                   - subject number of the current subject (necessary for filepaths and storage)
%   bemobil_config            - configuration struct with all necessary information. See EEG_processing_example file
%                                that comes with this function!
%
% Outputs:
%   ALLEEG                    - complete EEGLAB data set structure
%   EEG_preprocessed_and_ICA  - current EEGLAB EEG structure
%   Currentset                - index of current EEGLAB EEG structure within ALLEEG
%
%   .set data file of current EEGLAB EEG structure stored on disk 
%
% See also:
%   EEGLAB
%
% Authors: Marius Klug, 2021

function [ALLEEG, EEG_preprocessed_and_ICA, CURRENTSET] = bemobil_process_all_AMICA(ALLEEG, EEG_preprocessed, CURRENTSET,...
    subject, bemobil_config, force_recompute)

% check config
bemobil_config = bemobil_check_config(bemobil_config);

% get rid of memory mapped object storage and make sure double spacing and matlab save version 7 is used (for files
% larger than 2gb)
% mobilab uses memory mapped files which is why this needs to be set several times throughout the processing
try
    pop_editoptions( 'option_saveversion6', 0, 'option_single', 0, 'option_memmapdata', 0);
catch
    warning('Could NOT edit EEGLAB memory options!!');
end

if ~exist('force_recompute','var')
    force_recompute = false;
end
if force_recompute
    warning('RECOMPUTING OLD FILES IF FOUND!!!')
end

% check if the entire processing was done already
output_filepath = fullfile(bemobil_config.study_folder, bemobil_config.single_subject_analysis_folder,...
    [bemobil_config.filename_prefix num2str(subject)]);
mkdir(output_filepath)

if ~force_recompute
    try
        
        EEG = pop_loadset('filename', [bemobil_config.filename_prefix num2str(subject) '_'...
            bemobil_config.preprocessed_and_ICA_filename], 'filepath', output_filepath);
        [ALLEEG, EEG_preprocessed_and_ICA, CURRENTSET] = pop_newset(ALLEEG, EEG, 0,'study',0);
        
        warning('Old single subject file with interpolated channels, avref, and AMICA data already existed, skipping this processing!')
        
    catch
        disp('...failed. Computing now.')
    end
end

if ~exist('EEG_preprocessed_and_ICA','var')
    
    output_filepath = fullfile(bemobil_config.study_folder, bemobil_config.spatial_filters_folder,...
        bemobil_config.spatial_filters_folder_AMICA, [bemobil_config.filename_prefix num2str(subject)]);
    
    % check if the part was done already
    if ~force_recompute
        try
            
            EEG = pop_loadset('filename', [bemobil_config.filename_prefix num2str(subject) '_'...
                bemobil_config.amica_filename_output], 'filepath', output_filepath);
            [ALLEEG, EEG_AMICA, CURRENTSET] = pop_newset(ALLEEG, EEG, 0,'study',0);
            
            warning('Old AMICA file already existed, using that file!')
            
        catch
            disp('...failed. Computing now.')
        end
    end
    
    if ~exist('EEG_AMICA','var')
        
        %% highpass filter for AMICA
        
        [ALLEEG, EEG_filtered_for_AMICA, CURRENTSET] = bemobil_filter(ALLEEG, EEG_preprocessed, CURRENTSET,...
            bemobil_config.filter_lowCutoffFreqAMICA, bemobil_config.filter_highCutoffFreqAMICA,...
            [], [], bemobil_config.filter_AMICA_highPassOrder, bemobil_config.filter_AMICA_lowPassOrder);
        
        % save RAM and disk space of ICA results, since events are irrelevant here and in mobi datasets can be a lot
        EEG_filtered_for_AMICA.event = [];
        EEG_filtered_for_AMICA.urevent = [];
        
        %% AMICA
        % running signal decomposition with automatic rejection is recommended
        
        %         turns out the data rank reduction due to bridges does not really solve the issue of a high D value when
        %         computing AMICA, and since PCA before ICA was shown the be problematic, I removed it again: Artoni, F.,
        %         Delorme, A., & Makeig, S. (2018) Applying dimension reduction to EEG data by Principal Component Analysis
        %         reduces the quality of its subsequent Independent Component decomposition. Neuroimage, 175, 176-187.
        
        %         data_rank = EEG_filtered_for_AMICA.etc.rank;
        %         [rank_reduction_of_bridges,EEG_filtered_for_AMICA] = bemobil_find_gel_bridges(EEG_filtered_for_AMICA,0.98);
        %         data_rank = data_rank - rank_reduction_of_bridges;
        
        [ALLEEG, EEG_AMICA, CURRENTSET] = bemobil_signal_decomposition(ALLEEG, EEG_filtered_for_AMICA, ...
            CURRENTSET, true, bemobil_config.num_models, bemobil_config.max_threads, EEG_filtered_for_AMICA.etc.rank, [], ...
            [bemobil_config.filename_prefix num2str(subject) '_' bemobil_config.amica_filename_output], output_filepath,...
            bemobil_config.AMICA_autoreject, bemobil_config.AMICA_n_rej, bemobil_config.AMICA_reject_sigma_threshold,...
            bemobil_config.AMICA_max_iter);
        
        % plot autorejection
        data2plot = EEG_AMICA.data(1:round(EEG_AMICA.nbchan/10):EEG_AMICA.nbchan,:)';
        figure;
        set(gcf,'color','w','Position', get(0,'screensize'));
        plot(data2plot,'g');
        data2plot(~EEG_AMICA.etc.bad_samples,:) = NaN;
        hold on
        plot(data2plot,'r');
        xlim([-10000 EEG_AMICA.pnts+10000])
        ylim([-1000 1000])
        title(['AMICA autorejection, removed ' num2str(round(EEG_AMICA.etc.bad_samples_percent,2)) '% of the samples'])
        xlabel('Samples')
        ylabel('\muV')
        clear data2plot
        % save figure to disk
        savefig(gcf,fullfile(output_filepath,[bemobil_config.filename_prefix num2str(subject) '_AMICA_autoreject.fig']))
        print(gcf,fullfile(output_filepath,[bemobil_config.filename_prefix num2str(subject) '_AMICA_autoreject.png']),'-dpng')
        close
        
        % plot all ICs
        pop_topoplot(EEG_AMICA, 0, [1:size(EEG_AMICA.icaweights,1)],EEG_AMICA.filename,[],0,'electrodes','off');
        allICfighandle = gcf;
        print(allICfighandle,fullfile(output_filepath,[bemobil_config.filename_prefix num2str(subject) '_all_ICs.png']),'-dpng')
        close
    
        % save RAM
        clear EEG_filtered_for_AMICA
        
    end
    
    
    %% Warping of locations and dipole fitting, plus runing ICLabel
    % renames the specified channels, warps the chanlocs on a standard head model and fits dipoles for
    % each IC below the threshold of residual variance
    
    output_filepath = fullfile(bemobil_config.study_folder, bemobil_config.spatial_filters_folder,...
        bemobil_config.spatial_filters_folder_AMICA, [bemobil_config.filename_prefix num2str(subject)]);
    
    % check if the part was done already
    if ~force_recompute
        try
            
            EEG = pop_loadset('filename', [bemobil_config.filename_prefix num2str(subject) '_'...
                bemobil_config.dipfitted_filename], 'filepath', output_filepath);
            [ALLEEG, EEG_dipfitted, CURRENTSET] = pop_newset(ALLEEG, EEG, 0,'study',0);
            
            warning('Old dipfitted file already existed, using that file!')
            
        catch
            disp('...failed. Computing now.')
        end
    end
    
    if ~exist('EEG_dipfitted','var')
        
        % compute iclabel scores
        disp('ICLabel component classification...');
        EEG_AMICA = iclabel(EEG_AMICA,bemobil_config.iclabel_classifier);
        
        % do the warp and dipfit
        disp('Dipole fitting...');
        [ALLEEG, EEG_dipfitted, CURRENTSET] = bemobil_dipfit( EEG_AMICA , ALLEEG, CURRENTSET, bemobil_config.warping_channel_names,...
            bemobil_config.residualVariance_threshold,...
            bemobil_config.do_remove_outside_head, bemobil_config.number_of_dipoles,...
            [bemobil_config.filename_prefix num2str(subject) '_' bemobil_config.dipfitted_filename], output_filepath);
        
    end
    
    % save RAM
    clear EEG_AMICA
    
    %% Copy the spatial filter data into the raw full data set for further single subject processing
    
    output_filepath = fullfile(bemobil_config.study_folder, bemobil_config.single_subject_analysis_folder,...
        [bemobil_config.filename_prefix num2str(subject)]);
    
    disp('Copying all information into full length dataset for single subject processing...');
    [ALLEEG, EEG_single_subject_copied, CURRENTSET] = bemobil_copy_spatial_filter(EEG_preprocessed, ALLEEG, CURRENTSET,...
        EEG_dipfitted, [bemobil_config.filename_prefix num2str(subject) '_'...
        bemobil_config.preprocessed_and_ICA_filename], output_filepath);
    
    % save RAM
    clear EEG_preprocessed EEG_dipfitted
    %% clean with IClabel
    
    disp('Cleaning data with ICLabel')
    
    % clean now, save files and figs
    [ALLEEG, EEG_preprocessed_and_ICA, CURRENTSET, ICs_keep, ICs_throw] = bemobil_clean_with_iclabel( EEG_single_subject_copied ,...
        ALLEEG, CURRENTSET, bemobil_config.iclabel_classifier,...
        bemobil_config.iclabel_classes, bemobil_config.iclabel_threshold,...
        [ bemobil_config.filename_prefix num2str(subject) '_' bemobil_config.single_subject_cleaned_ICA_filename],output_filepath);
    
    disp('...done.')
    
    %% plot cleaned with ICA, for analytics
    
    disp('Filtering data only for plotting!')
    EEG = pop_eegfiltnew(EEG_preprocessed_and_ICA, 'locutoff',0.5);
    
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
    axeshandle = gca;
    fighandle = gcf;
    axcp = copyobj(axeshandle, plotfigure);
    set(axcp,'Position',get(ax1,'position'));
    axcp.XTickLabel = [0:10]+round(starttime/1000);
    axcp.YTick=[];
    axcp.Title.String = ['Cleaned channels data section 1 of ' num2str(round(EEG.times(end)/1000)) 's'];
    axcp.XLabel.String = 'seconds';
    delete(ax1);
    close(fighandle)
    
    starttime = EEG.times(end)/7*2;
    vis_artifacts(EEG,EEG,'show_events',1,'time_subset',...
        round([starttime starttime+10000]/1000)); % plot 10s at the first quarter
    axeshandle = gca;
    fighandle = gcf;
    axcp = copyobj(axeshandle, plotfigure);
    set(axcp,'Position',get(ax2,'position'));
    axcp.XTickLabel = [0:10]+round(starttime/1000);
    axcp.YTick=[];
    axcp.Title.String = ['Cleaned channels data section 2 of ' num2str(round(EEG.times(end)/1000)) 's'];
    axcp.XLabel.String = 'seconds';
    delete(ax2);
    close(fighandle)
    
    starttime = EEG.times(end)/7*3;
    vis_artifacts(EEG,EEG,'show_events',1,'time_subset',...
        round([starttime starttime+10000]/1000)); % plot 10s at the first quarter
    axeshandle = gca;
    fighandle = gcf;
    axcp = copyobj(axeshandle, plotfigure);
    set(axcp,'Position',get(ax3,'position'));
    axcp.XTickLabel = [0:10]+round(starttime/1000);
    axcp.YTick=[];
    axcp.Title.String = ['Cleaned channels data section 3 of ' num2str(round(EEG.times(end)/1000)) 's'];
    axcp.XLabel.String = 'seconds';
    delete(ax3);
    close(fighandle)
    
    starttime = EEG.times(end)/7*4;
    vis_artifacts(EEG,EEG,'show_events',1,'time_subset',...
        round([starttime starttime+10000]/1000)); % plot 10s at the first quarter
    axeshandle = gca;
    fighandle = gcf;
    axcp = copyobj(axeshandle, plotfigure);
    set(axcp,'Position',get(ax4,'position'));
    axcp.XTickLabel = [0:10]+round(starttime/1000);
    axcp.YTick=[];
    axcp.Title.String = ['Cleaned channels data section 4 of ' num2str(round(EEG.times(end)/1000)) 's'];
    axcp.XLabel.String = 'seconds';
    delete(ax4);
    close(fighandle)
    
    starttime = EEG.times(end)/7*5;
    vis_artifacts(EEG,EEG,'show_events',1,'time_subset',...
        round([starttime starttime+10000]/1000)); % plot 10s at the first quarter
    axeshandle = gca;
    fighandle = gcf;
    axcp = copyobj(axeshandle, plotfigure);
    set(axcp,'Position',get(ax5,'position'));
    axcp.XTickLabel = [0:10]+round(starttime/1000);
    axcp.YTick=[];
    axcp.Title.String = ['Cleaned channels data section 5 of ' num2str(round(EEG.times(end)/1000)) 's'];
    axcp.XLabel.String = 'seconds';
    delete(ax5);
    close(fighandle)
    
    starttime = EEG.times(end)/7*6;
    vis_artifacts(EEG,EEG,'show_events',1,'time_subset',...
        round([starttime starttime+10000]/1000)); % plot 10s at the first quarter
    axeshandle = gca;
    fighandle = gcf;
    axcp = copyobj(axeshandle, plotfigure);
    set(axcp,'Position',get(ax6,'position'));
    axcp.XTickLabel = [0:10]+round(starttime/1000);
    axcp.YTick=[];
    axcp.Title.String = ['Cleaned channels data section 6 of ' num2str(round(EEG.times(end)/1000)) 's'];
    axcp.XLabel.String = 'seconds';
    delete(ax6);
    close(fighandle)
    
    %% save plot
    
    savefig(plotfigure,fullfile(output_filepath,[bemobil_config.filename_prefix num2str(subject) '_cleaned.fig']))
    print(plotfigure,fullfile(output_filepath,[bemobil_config.filename_prefix num2str(subject) '_cleaned.png']),'-dpng')
    close
    
end

disp('Entire AMICA processing done!');
