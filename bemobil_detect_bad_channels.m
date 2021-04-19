% bemobil_detect_bad_channels - Detects bad channels using the "clean_artifacts" function of the EEGLAB "clean_rawdata"
% plugin. Computes average reference before detecting bad channels and uses an 0.5Hz highpass filter cutoff (inside
% clan_artifacts). Plots segments of the data to check the detection.
%
% Usage:
%   >>  [chans_to_interp, plothandle] = bemobil_detect_bad_channels(EEG, ALLEEG, CURRENTSET, chancorr_crit, chan_max_broken_time)
% 
% Inputs:
%   EEG                     - current EEGLAB EEG structure
%   ALLEEG                  - complete EEGLAB data set structure
%   CURRENTSET              - index of current EEGLAB EEG structure within ALLEEG
%   chancorr_crit           - Correlation threshold. If a channel is correlated at less than this value
%                               to its robust estimate (based on other channels), it is considered abnormal in
%                               the given time window. OPTIONAL, default = 0.8.
%   chan_max_broken_time    - Maximum time (either in seconds or as fraction of the recording) during which a 
%                               retained channel may be broken. Reasonable range: 0.1 (very aggressive) to 0.6
%                               (very lax). OPTIONAL, default = 0.5.
%
% Outputs:
%   chans_to_interp         - vector with channel indices to remove
%   plothandle              - handle to the plot of the data segments to check cleaning
%
%   .set data file of current EEGLAB EEG structure stored on disk (OPTIONALLY)
%
% See also:
%   EEGLAB, bemobil_avref, clean_artifacts
%
% Authors: Lukas Gehrke, 2017, Marius Klug, 2021

function [chans_to_interp, plothandle] = bemobil_detect_bad_channels(EEG, ALLEEG, CURRENTSET, chancorr_crit, chan_max_broken_time)

if ~exist('chancorr_crit','var') || isempty(chancorr_crit)
	chancorr_crit = 0.8;
end
if ~exist('chan_max_broken_time','var') || isempty(chan_max_broken_time)
	chan_max_broken_time = 0.5;
end

if ~strcmp(EEG.ref,'average')
    disp('Re-referencing ONLY for bad channel detection now.')
    % compute average reference before finding bad channels 
    [ALLEEG, EEG, CURRENTSET] = bemobil_avref( EEG , ALLEEG, CURRENTSET);
end

disp('Detecting bad channels...')

% remove bad channels, use default values of clean_artifacts, but specify just in case they may change
[EEG_chan_removed,EEG_highpass,~,chans_to_interp] = clean_artifacts(EEG,...
    'burst_crit','off','window_crit','off','ChannelCriterionMaxBadTime',chan_max_broken_time,...
    'chancorr_crit',chancorr_crit,'line_crit',4,'highpass_band',[0.25 0.75],'flatline_crit','on');
chans_to_interp = find(chans_to_interp); % transform logical array to indices

disp('Detected bad channels: ')
disp({EEG.chanlocs(chans_to_interp).labels})

% take EOG out of the channels to interpolate: EOG is very likely to be different from the others but rightfully so
disp('Ignoring EOG channels for interpolation:')

disp({EEG.chanlocs(chans_to_interp(strcmp({EEG.chanlocs(chans_to_interp).type},'EOG'))).labels})
chans_to_interp(strcmp({EEG.chanlocs(chans_to_interp).type},'EOG'))=[];

disp('Final bad channels: ')
disp({EEG.chanlocs(chans_to_interp).labels})

EEG_chan_removed = pop_select( EEG_highpass,'nochannel',chans_to_interp);
EEG_chan_removed.etc.clean_channel_mask = true(EEG_highpass.nbchan,1);
EEG_chan_removed.etc.clean_channel_mask(chans_to_interp) = deal(0);

% display 1/10 of the data in the middle (save disk space when saving figure)


plothandle = figure('color','w');
set(plothandle, 'Position', get(0,'screensize'))
ax1 = subplot(231);
ax2 = subplot(232);
ax3 = subplot(233);
ax4 = subplot(234);
ax5 = subplot(235);
ax6 = subplot(236);


starttime = EEG.times(end)/7*1;
vis_artifacts(EEG_chan_removed,EEG,'show_events',1,'time_subset',...
    round([starttime starttime+10000]/1000)); % plot 10s at the first quarter
axeshandle = gca;
fighandle = gcf;
axcp = copyobj(axeshandle, plothandle);
set(axcp,'Position',get(ax1,'position'));
axcp.XTickLabel = [0:10]+round(starttime/1000);
axcp.YTick=[];
axcp.Title.String = ['Bad channels data section 1 of ' num2str(round(EEG.times(end)/1000)) 's'];
axcp.XLabel.String = 'seconds';
delete(ax1);
close(fighandle)

starttime = EEG.times(end)/7*2;
vis_artifacts(EEG_chan_removed,EEG,'show_events',1,'time_subset',...
    round([starttime starttime+10000]/1000)); % plot 10s at the first quarter
axeshandle = gca;
fighandle = gcf;
axcp = copyobj(axeshandle, plothandle);
set(axcp,'Position',get(ax2,'position'));
axcp.XTickLabel = [0:10]+round(starttime/1000);
axcp.YTick=[];
axcp.Title.String = ['Bad channels data section 2 of ' num2str(round(EEG.times(end)/1000)) 's'];
axcp.XLabel.String = 'seconds';
delete(ax2);
close(fighandle)

starttime = EEG.times(end)/7*3;
vis_artifacts(EEG_chan_removed,EEG,'show_events',1,'time_subset',...
    round([starttime starttime+10000]/1000)); % plot 10s at the first quarter
axeshandle = gca;
fighandle = gcf;
axcp = copyobj(axeshandle, plothandle);
set(axcp,'Position',get(ax3,'position'));
axcp.XTickLabel = [0:10]+round(starttime/1000);
axcp.YTick=[];
axcp.Title.String = ['Bad channels data section 3 of ' num2str(round(EEG.times(end)/1000)) 's'];
axcp.XLabel.String = 'seconds';
delete(ax3);
close(fighandle)

starttime = EEG.times(end)/7*4;
vis_artifacts(EEG_chan_removed,EEG,'show_events',1,'time_subset',...
    round([starttime starttime+10000]/1000)); % plot 10s at the first quarter
axeshandle = gca;
fighandle = gcf;
axcp = copyobj(axeshandle, plothandle);
set(axcp,'Position',get(ax4,'position'));
axcp.XTickLabel = [0:10]+round(starttime/1000);
axcp.YTick=[];
axcp.Title.String = ['Bad channels data section 4 of ' num2str(round(EEG.times(end)/1000)) 's'];
axcp.XLabel.String = 'seconds';
delete(ax4);
close(fighandle)

starttime = EEG.times(end)/7*5;
vis_artifacts(EEG_chan_removed,EEG,'show_events',1,'time_subset',...
    round([starttime starttime+10000]/1000)); % plot 10s at the first quarter
axeshandle = gca;
fighandle = gcf;
axcp = copyobj(axeshandle, plothandle);
set(axcp,'Position',get(ax5,'position'));
axcp.XTickLabel = [0:10]+round(starttime/1000);
axcp.YTick=[];
axcp.Title.String = ['Bad channels data section 5 of ' num2str(round(EEG.times(end)/1000)) 's'];
axcp.XLabel.String = 'seconds';
delete(ax5);
close(fighandle)

starttime = EEG.times(end)/7*6;
vis_artifacts(EEG_chan_removed,EEG,'show_events',1,'time_subset',...
    round([starttime starttime+10000]/1000)); % plot 10s at the first quarter
axeshandle = gca;
fighandle = gcf;
axcp = copyobj(axeshandle, plothandle);
set(axcp,'Position',get(ax6,'position'));
axcp.XTickLabel = [0:10]+round(starttime/1000);
axcp.YTick=[];
axcp.Title.String = ['Bad channels data section 6 of ' num2str(round(EEG.times(end)/1000)) 's'];
axcp.XLabel.String = 'seconds';
delete(ax6);
close(fighandle)