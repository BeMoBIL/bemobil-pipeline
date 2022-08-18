% bemobil_detect_bad_channels - Repeatedly detects bad channels using the "clean_artifacts" function of the EEGLAB
% "clean_rawdata" plugin. Computes average reference before detecting bad channels and uses an 0.5Hz highpass filter
% cutoff (inside clan_artifacts). Plots the rejections of each iteration and the final rejections. Plots segments of the
% data to check the rejected channels.
%
% Usage:
%   >>  [chans_to_interp, chan_detected_fraction_threshold, detected_bad_channels, rejected_chan_plot_handle, detection_plot_handle] =...
%           bemobil_detect_bad_channels(EEG, ALLEEG, CURRENTSET, chancorr_crit, chan_max_broken_time, chan_detect_num_iter,...
%           chan_detected_fraction_threshold, num_chan_rej_max_target, flatline_crit, line_noise_crit)
% 
% Inputs:
%   EEG                                 - current EEGLAB EEG structure
%   ALLEEG                              - complete EEGLAB data set structure
%   CURRENTSET                          - index of current EEGLAB EEG structure within ALLEEG
%   chancorr_crit                       - Correlation threshold. If a channel is correlated at less than this value
%                                           to its robust estimate (based on other channels), it is considered abnormal in
%                                           the given time window. OPTIONAL, default = 0.8.
%   chan_max_broken_time                - Maximum time (either in seconds or as fraction of the recording) during which a 
%                                           retained channel may be broken. Reasonable range: 0.1 (very aggressive) to 0.6
%                                           (very lax). OPTIONAL, default = 0.5.
%   chan_detect_num_iter                - Number of iterations the bad channel detection should run (default = 10)
%   chan_detected_fraction_threshold	- Fraction how often a channel has to be detected to be rejected in the final
%                                           rejection (default 0.5)
%   num_chan_rej_max_target             - Target max amount of channel rejection. Actual num of rejections might be higher if 
%                                           there are a lot of bad channels. Target precision can be increased with higher chan_detect_num_iter
%                                           If empty, use only chan_detected_fraction_threshold. Can be either a fraction 
%                                           of all channels (will be rounded, e.g. 1/5 of chans) or a specific integer
%                                           number. (default = 1/5)
%   flatline_crit                       - Maximum duration a channel can be flat in seconds (default 'off')
%   line_noise_crit                     - If a channel has more line noise relative to its signal than this value, in
%                                           standard deviations based on the total channel population, it is considered
%                                           abnormal. (default: 'off')
%
% Outputs:
%   chans_to_interp                     - vector with channel indices to remove
%   chan_detected_fraction_threshold    - 'badness' threshold used for bad channel selection. Might deviate from 
%                                           chan_detected_fraction_threshold when there are a lot of bad channels.
%   detected_bad_channels               - n*m matrix of boolean indicating what channel was marked for rejection in which iteration
%                                           n is number of channels and m number of iterations
%   rejected_chan_plot_handle           - handle to the plot of the data segments to check cleaning
%   detection_plot_handle               - handle to the plot of the rejection iterations
%   
%
%   .set data file of current EEGLAB EEG structure stored on disk (OPTIONALLY)
%
% See also:
%   EEGLAB, bemobil_avref, clean_artifacts
%
% Authors: Lukas Gehrke, 2017, Marius Klug, 2021, Timotheus Berg, 2021

function [chans_to_interp, chan_detected_fraction_threshold, detected_bad_channels, rejected_chan_plot_handle, detection_plot_handle] = bemobil_detect_bad_channels(EEG, ALLEEG, CURRENTSET,...
    chancorr_crit, chan_max_broken_time, chan_detect_num_iter, chan_detected_fraction_threshold, num_chan_rej_max_target, flatline_crit, line_noise_crit)

if ~exist('chancorr_crit','var') || isempty(chancorr_crit)
	chancorr_crit = 0.8;
end
if ~exist('chan_max_broken_time','var') || isempty(chan_max_broken_time)
	chan_max_broken_time = 0.5;
end
if ~exist('chan_detect_num_iter','var') || isempty(chan_detect_num_iter)
	chan_detect_num_iter = 10;
end
if ~exist('chan_detected_fraction_threshold','var') || isempty(chan_detected_fraction_threshold)
	chan_detected_fraction_threshold = 0.5;
end
if ~exist('flatline_crit','var') || isempty(flatline_crit)
	flatline_crit = 'off';
end
if ~exist('line_noise_crit','var') || isempty(line_noise_crit)
	line_noise_crit = 'off';
end
if ~exist('num_chan_rej_target','var') || isempty(num_chan_rej_max_target)
	num_chan_rej_max_target = 1/5;
end

%%
if ~strcmp(EEG.ref,'average')
    disp('Re-referencing ONLY for bad channel detection now.')
    % compute average reference before finding bad channels 
    [ALLEEG, EEG, CURRENTSET] = bemobil_avref( EEG , ALLEEG, CURRENTSET);
end

%%
disp('Repeated bad channels detection...')
detected_bad_channels = [];
for i = 1:chan_detect_num_iter
    disp(['Iteration ' num2str(i) '/' num2str(chan_detect_num_iter)])
    clear hlp_microcache
    % remove bad channels, use default values of clean_artifacts, but specify just in case they may change
    [EEG_chan_removed,EEG_highpass,~,detected_bad_channels(1:EEG.nbchan,i)] = clean_artifacts(EEG,...
        'burst_crit','off','window_crit','off','channel_crit_maxbad_time',chan_max_broken_time,...
        'chancorr_crit',chancorr_crit,'line_crit',line_noise_crit,'highpass_band',[0.25 0.75],'flatline_crit',flatline_crit);

end
disp('...iterative bad channel detection done!')

badness_percent = sum(detected_bad_channels,2) / size(detected_bad_channels,2);

% If there are a lot of bad channels, raise the threshold to the badness value at num_chan_rej_target
if ~isempty(num_chan_rej_max_target)
    if num_chan_rej_max_target<1 && num_chan_rej_max_target>0
        num_chan_rej_max_target = round(EEG.nbchan * num_chan_rej_max_target);
    end
    [sort_val, ~ ] = sort(badness_percent, 'descend');

    if sort_val(num_chan_rej_max_target) > chan_detected_fraction_threshold
        chan_detected_fraction_threshold = sort_val(num_chan_rej_max_target);
    end
end

chans_to_interp = badness_percent >= chan_detected_fraction_threshold;

%% plot
detection_plot_handle = figure; 
set(detection_plot_handle,'color','w','position',[501 325 1352 741])

subplot(1,50,[1:39])
imagesc(detected_bad_channels);
colormap cool
title(['detected bad channels, necessary fraction for removal: ' num2str(chan_detected_fraction_threshold)])
xlabel('iteration')
ylabel('channel')
set(gca,'fontsize',12)

subplot(1,50,[42:46])
imagesc(badness_percent);
xticks([])
yticks([])
title('badness')
set(gca,'fontsize',12)
colorbar

subplot(1,50,[48:50])
imagesc(chans_to_interp);
xticks([])
yticks([])
title('final')
set(gca,'fontsize',12)

drawnow

%% select the final channels to remove and remove them from a dataset to plot

% give actual channel numbers as output
chans_to_interp = find(chans_to_interp);

disp('Detected bad channels: ')
disp({EEG.chanlocs(chans_to_interp).labels})

% take EOG out of the channels to interpolate: EOG is very likely to be different from the others but rightfully so
disp('Ignoring EOG channels for interpolation:')

disp({EEG.chanlocs(chans_to_interp(strcmp({EEG.chanlocs(chans_to_interp).type},'EOG'))).labels})
chans_to_interp(strcmp({EEG.chanlocs(chans_to_interp).type},'EOG'))=[];

disp('Final bad channels: ')
disp({EEG.chanlocs(chans_to_interp).labels})


% remove channels and store channel mask for plotting
EEG_chan_removed = pop_select( EEG_highpass,'nochannel',chans_to_interp);

clean_channel_mask = ones(EEG.nbchan,1);
clean_channel_mask(chans_to_interp) = 0;
EEG_chan_removed.etc.clean_channel_mask = logical(clean_channel_mask);

%% plot

rejected_chan_plot_handle = figure('color','w');
set(rejected_chan_plot_handle, 'Position', get(0,'screensize'))
ax1 = subplot(231);
ax2 = subplot(232);
ax3 = subplot(233);
ax4 = subplot(234);
ax5 = subplot(235);
ax6 = subplot(236);

titletext = [];
for i = 1:length(chans_to_interp)
    titletext = [titletext EEG.chanlocs(chans_to_interp(i)).labels ' (' num2str(chans_to_interp(i)) '), '];
end
titletext = titletext(1:end-2);

h = textsc(['Detected: ' titletext], 'title');
set(h, 'fontsize',16)

starttime = EEG.times(end)/7*1;
vis_artifacts(EEG_chan_removed,EEG,'show_events',1,'time_subset',...
    round([starttime starttime+10000]/1000),'equalize_channel_scaling',1); % plot 10s at the first quarter
drawnow
axeshandle = gca;
fighandle = gcf;
axcp = copyobj(axeshandle, rejected_chan_plot_handle);
set(axcp,'Position',get(ax1,'position'));
axcp.XTickLabel = [0:10]+round(starttime/1000);
axcp.YTick=[];
axcp.Title.String = ['Bad channels data section 1 of ' num2str(round(EEG.times(end)/1000)) 's'];
axcp.XLabel.String = 'seconds';
drawnow
delete(ax1);
close(fighandle)

starttime = EEG.times(end)/7*2;
vis_artifacts(EEG_chan_removed,EEG,'show_events',1,'time_subset',...
    round([starttime starttime+10000]/1000),'equalize_channel_scaling',1); % plot 10s at the first quarter
drawnow
axeshandle = gca;
fighandle = gcf;
axcp = copyobj(axeshandle, rejected_chan_plot_handle);
set(axcp,'Position',get(ax2,'position'));
axcp.XTickLabel = [0:10]+round(starttime/1000);
axcp.YTick=[];
axcp.Title.String = ['Bad channels data section 2 of ' num2str(round(EEG.times(end)/1000)) 's'];
axcp.XLabel.String = 'seconds';
drawnow
delete(ax2);
close(fighandle)

starttime = EEG.times(end)/7*3;
vis_artifacts(EEG_chan_removed,EEG,'show_events',1,'time_subset',...
    round([starttime starttime+10000]/1000),'equalize_channel_scaling',1); % plot 10s at the first quarter
drawnow
axeshandle = gca;
fighandle = gcf;
axcp = copyobj(axeshandle, rejected_chan_plot_handle);
set(axcp,'Position',get(ax3,'position'));
axcp.XTickLabel = [0:10]+round(starttime/1000);
axcp.YTick=[];
axcp.Title.String = ['Bad channels data section 3 of ' num2str(round(EEG.times(end)/1000)) 's'];
axcp.XLabel.String = 'seconds';
drawnow
delete(ax3);
close(fighandle)

starttime = EEG.times(end)/7*4;
vis_artifacts(EEG_chan_removed,EEG,'show_events',1,'time_subset',...
    round([starttime starttime+10000]/1000),'equalize_channel_scaling',1); % plot 10s at the first quarter
drawnow
axeshandle = gca;
fighandle = gcf;
axcp = copyobj(axeshandle, rejected_chan_plot_handle);
set(axcp,'Position',get(ax4,'position'));
axcp.XTickLabel = [0:10]+round(starttime/1000);
axcp.YTick=[];
axcp.Title.String = ['Bad channels data section 4 of ' num2str(round(EEG.times(end)/1000)) 's'];
axcp.XLabel.String = 'seconds';
drawnow
delete(ax4);
close(fighandle)

starttime = EEG.times(end)/7*5;
vis_artifacts(EEG_chan_removed,EEG,'show_events',1,'time_subset',...
    round([starttime starttime+10000]/1000),'equalize_channel_scaling',1); % plot 10s at the first quarter
drawnow
axeshandle = gca;
fighandle = gcf;
axcp = copyobj(axeshandle, rejected_chan_plot_handle);
set(axcp,'Position',get(ax5,'position'));
axcp.XTickLabel = [0:10]+round(starttime/1000);
axcp.YTick=[];
axcp.Title.String = ['Bad channels data section 5 of ' num2str(round(EEG.times(end)/1000)) 's'];
axcp.XLabel.String = 'seconds';
drawnow
delete(ax5);
close(fighandle)

starttime = EEG.times(end)/7*6;
vis_artifacts(EEG_chan_removed,EEG,'show_events',1,'time_subset',...
    round([starttime starttime+10000]/1000),'equalize_channel_scaling',1); % plot 10s at the first quarter
drawnow
axeshandle = gca;
fighandle = gcf;
axcp = copyobj(axeshandle, rejected_chan_plot_handle);
set(axcp,'Position',get(ax6,'position'));
axcp.XTickLabel = [0:10]+round(starttime/1000);
axcp.YTick=[];
axcp.Title.String = ['Bad channels data section 6 of ' num2str(round(EEG.times(end)/1000)) 's'];
axcp.XLabel.String = 'seconds';
drawnow
delete(ax6);
close(fighandle)

% in case one clicks on the figure while processing the colormap gets lost so here it is again
figure(detection_plot_handle)
colormap cool