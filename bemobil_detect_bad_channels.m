function [chans_to_interp, plothandle] = bemobil_detect_bad_channels(EEG, ALLEEG, CURRENTSET, chancorr_crit)

if ~strcmp(EEG.ref,'average')
    % compute average reference before finding bad channels 
    [ALLEEG, EEG, CURRENTSET] = bemobil_avref( EEG , ALLEEG, CURRENTSET);
end

disp('Detecting bad channels...')

% remove bad channels, use default values of clean_artifacts, but specify just in case they may change
[EEG_chan_removed,EEG_highpass,~,chans_to_interp] = clean_artifacts(EEG,...
    'burst_crit','off','window_crit','off',...
    'chancorr_crit',chancorr_crit,'line_crit',4,'highpass_band',[0.25 0.75],'flatline_crit','off');
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
vis_artifacts(EEG_chan_removed,EEG,'show_events',0,'time_subset',...
    [round(EEG.times(end)/2) round(EEG.times(end)/2+round(EEG.times(end)/10))]/1000);

plothandle = gcf;
set(plothandle, 'Position', get(0,'screensize'))