% bemobil_preprocess() - Preprocessing of EEG data: Fill EEG structure with
% ur-data. Remove unused electrodes of electrode arrays. Import channel
% locations from vicra file. Change channel names and declare EOG.
% Resample.
%
% Usage:
%   >>  [ ALLEEG EEG CURRENTSET ] = bemobil_preprocess(ALLEEG, EEG, CURRENTSET, channel_locations_filepath, channels_to_remove, eog_channels, resample_freq);
%   >>  [ ALLEEG EEG CURRENTSET ] = bemobil_preprocess(ALLEEG, EEG, CURRENTSET, channel_locations_filepath, channels_to_remove, eog_channels, resample_freq, out_filename, out_filepath);
%
% Inputs:
%   ALLEEG                  - complete EEGLAB data set structure
%   EEG                     - current EEGLAB EEG structure
%   CURRENTSET              - index of current EEGLAB EEG structure within ALLEEG
%   channel_locations_file  - channel_locations file (with path); OR []
%   channels_to_remove      - cell of all channels that should be thrown out
%       per se (e.g. {'N29' 'N30' 'N31'}); OR []
%   eog_channels            - cell of channels that should be declared as EOG
%       for later use (e.g. {'G16' 'G32'}); OR []
%   resample_freq           - Resample frequency (Hz), if [], no resampling will be applied
%   out_filename            - output filename (OPTIONAL ARGUMENT)
%   out_filepath            - output filepath (OPTIONAL ARGUMENT - File will only be saved on disk
%       if both a name and a path are provided)
%
% Outputs:
%   ALLEEG                  - complete EEGLAB data set structure
%   EEG                     - current EEGLAB EEG structure
%   Currentset              - index of current EEGLAB EEG structure within ALLEEG
%
%   .set data file of current EEGLAB EEG structure stored on disk (OPTIONALLY)
%
% See also:
%   EEGLAB, pop_eegfiltnew, pop_resample, pop_chanedit, pop_select
%
% Authors: Lukas Gehrke, 2017, Marius Klug, 2021

function [ ALLEEG EEG CURRENTSET ] = bemobil_preprocess(ALLEEG, EEG, CURRENTSET,...
    channel_locations_filepath, channels_to_remove, eog_channels, resample_freq,...
    out_filename, out_filepath, rename_channels, ref_channel, linefreqs, zapline_n_remove, plot_zapline)

% only save a file on disk if both a name and a path are provided
save_file_on_disk = (exist('out_filename', 'var') && ~isempty(out_filename) &&...
    exist('out_filepath', 'var') && ~isempty(out_filepath));

% check if file already exist and show warning if it does
if save_file_on_disk
    mkdir(out_filepath); % make sure that folder exists, nothing happens if so
    dir_files = dir(out_filepath);
    if ismember(out_filename, {dir_files.name})
        warning([out_filename ' file already exists in: ' out_filepath '. File will be overwritten...']);
    end
end

% fill/copy all ur_structures with raw data (e.g. copy event to urevent)
EEG = eeg_checkset(EEG, 'makeur');

% remove unused neck electrodes from file (if BeMoBIL layout is used as is)
if ~isempty(channels_to_remove)
    if all(ismember(channels_to_remove, {EEG.chanlocs.labels}))
        % b) remove not needed channels import "corrected" chanlocs file
        EEG = pop_select( EEG,'nochannel', channels_to_remove);
        EEG = eeg_checkset( EEG );
        disp(['Removed electrodes: ' channels_to_remove ' from the dataset.']);
    else
        error('Not all of the specified channels to remove were present als data channels!')
    end
else
    disp('No channels to remove specified, skipping this step.')
end

% Resample/downsample to 250 Hz if no other resampling frequency is
% provided
if ~isempty(resample_freq)
    EEG = pop_resample(EEG, resample_freq);
    EEG = eeg_checkset( EEG );
    disp(['Resampled data to: ', num2str(resample_freq), 'Hz.']);
end

%% Clean line noise with ZapLine: de Cheveigné, A. (2020) ZapLine: a simple and effective method to remove power line
% artifacts. Neuroimage, 1, 1–13.
if exist('linefreqs','var') && ~isempty(linefreqs)
    
    disp('Removing frequency artifacts using ZapLine with adaptations for automatic component selection.')
    disp('---------------- PLEASE CITE ------------------')
    disp('de Cheveigné, A. (2020) ZapLine: a simple and effective method to remove power line artifacts. Neuroimage, 1, 1–13.')
    disp('---------------- PLEASE CITE ------------------')
    
    if ~exist('zapline_n_remove','var') || isempty(zapline_n_remove)
        disp('No nremove value was defined, using adaptive threshold detection for ZapLine (recommended).')
        zapline_config.adaptiveNremove = 1;
        zapline_n_remove_start = 1;
    else
        fprintf('ZapLine nremove value was defined, removing %d components.\n', zapline_n_remove);
        warning('Note: using adaptive threshold detection for ZapLine is recommended. For this, remove the zapline_n_remove parameter.')
        zapline_config.adaptiveNremove = 0;
        zapline_n_remove_start = zapline_n_remove;
    end
    
     if ~exist('plot_zapline','var') || isempty(plot_zapline)
         plot_zapline = 1;
     end
    
    % Step through each frequency one after another
    for i_linefreq = 1:length(linefreqs)
        %         error('nope')
        linefreq = linefreqs(i_linefreq);
        fline = linefreq/EEG.srate;
        
        fprintf('Removing noise at %gHz... \n',linefreq);
        
        zapline_config.fig2 = 100+i_linefreq;
        zapline_config.nfft = 512;
        %         zapline_config.nkeep = 40;
        
        if 0
            % full dataset zapline
            
            p.fig1 = 10;
            p.fig2 = 11;
            p.nfft = 512;
            
            [y,~,zapline_n_remove_final, scores] = nt_zapline_bemobil(EEG.data',fline,zapline_n_remove_start,zapline_config,0);
            
            x = EEG.data';
            
            figure(p.fig1);clf;set(gcf,'color','w')
            plot(scores);
            hold on
            ylim([0,max(scores)])
            xlim([0, length(scores)])
            plot([zapline_n_remove_final zapline_n_remove_final],ylim)
            title(['Artifact scores of ' num2str(linefreq) 'Hz components, removed = ' num2str(zapline_n_remove_final)])
            xlabel('Component')
            drawnow
            
            disp('Done. Computing spectra for plotting.')
            
            figure(p.fig2)
            clf; set(gcf,'color','w','Position',[1 1 1920 1080])
            figure(p.fig2)
            subplot 121
            [pxx,f]=nt_spect_plot(x/sqrt(mean(x(:).^2)),p.nfft,[],[],1/fline);
            divisor=sum(pxx);
            semilogy(f,abs(pxx)/divisor);
            legend('original'); legend boxoff
            set(gca,'ygrid','on','xgrid','on');
            xlabel('frequency (relative to line)');
            ylabel('relative power');
            yl1=get(gca,'ylim');
            hh=get(gca,'children');
            set(hh(1),'color','k')
            drawnow
            
            figure(p.fig2)
            subplot 122
            [pxx,f]=nt_spect_plot(y/sqrt(mean(x(:).^2)),p.nfft,[],[],1/fline);
            figure(p.fig2)
            semilogy(f,abs(pxx)/divisor);
            drawnow
            
            hold on
            [pxx,f]=nt_spect_plot((x-y)/sqrt(mean(x(:).^2)),p.nfft,[],[],1/fline);
            figure(p.fig2)
            semilogy(f,abs(pxx)/divisor);
            legend('clean', 'removed'); legend boxoff
            set(gca,'ygrid','on','xgrid','on');
            set(gca,'yticklabel',[]); ylabel([]);
            xlabel('frequency (relative to line)');
            yl2=get(gca,'ylim');
            hh=get(gca,'children');
            set(hh(1),'color',[1 .5 .5]); set(hh(2), 'color', [ 0 .7 0]);
            set(hh(2),'linewidth', 2);
            yl(1)=min(yl1(1),yl2(1)); yl(2)=max(yl1(2),yl2(2));
            subplot 121; ylim(yl); subplot 122; ylim(yl);
            drawnow
            
            disp('...done.')
            
            
            
        else
            % zapline on chunked data
            
            x = EEG.data';
            y = NaN(size(x));
            length_chunks = 60;
            n_chunks = floor(size(x,1)/250/length_chunks); % last chunk must be larger than the others, to ensure fft works
            scores = NaN(n_chunks,size(x,2));
            for i_chunk = 1:n_chunks
                if mod(i_chunk,round(n_chunks/10))==0
                    disp(['Chunk ' num2str(i_chunk) ' of ' num2str(n_chunks)])
                end
                
                if i_chunk ~= n_chunks
                    chunk_indices = 1+length_chunks*EEG.srate*(i_chunk-1):length_chunks*EEG.srate*(i_chunk);
                else
                    chunk_indices = 1+length_chunks*EEG.srate*(i_chunk-1):EEG.pnts;
                end
                
                chunk = x(chunk_indices,:);
                [y(chunk_indices,:),...
                    ~,zapline_n_remove_final(i_chunk), this_scores] =...
                    nt_zapline_bemobil(chunk,fline,zapline_n_remove_start,zapline_config,0);
                scores(i_chunk,1:length(this_scores)) = this_scores;
            end
            
            if plot_zapline
                disp('Done. Computing spectra for plotting.')
                
                % plot
                
                p.fig1 = 100;
                p.fig2 = 101;
                p.nfft = 512;
                
                figure(p.fig1);
                clf; set(gcf,'color','w','Position',[1 1 1920 420])
                subplot(1,5,[1:4]);
                plot(zapline_n_remove_final)
                xlabel('Chunk')
                ylabel('Removed Comps')
                title(['Removed ' num2str(linefreq) 'Hz components per ' num2str(length_chunks) 's chunk (automatic detection), M = ' num2str(mean(zapline_n_remove_final))])
                subplot(155)
                plot(nanmean(scores))
                title(['Mean artifact scores'])
                xlabel('Component')
                
                disp('proportion of non-DC power removed:');
                disp(nt_wpwr(x-y)/nt_wpwr(nt_demean(x)));
                
                figure(p.fig2);
                clf; set(gcf,'color','w','Position',[1 1 1920 1080])
                subplot 121
                [pxx,f]=nt_spect_plot(x/sqrt(mean(x(:).^2)),p.nfft,[],[],1/fline);
                divisor=sum(pxx);
                figure(p.fig2);
                semilogy(f,abs(pxx)/divisor);
                legend('original'); legend boxoff
                set(gca,'ygrid','on','xgrid','on');
                xlabel('frequency (relative to line)');
                ylabel('relative power');
                yl1=get(gca,'ylim');
                hh=get(gca,'children');
                set(hh(1),'color','k')
                drawnow
                subplot 122
                [pxx,f]=nt_spect_plot(y/sqrt(mean(x(:).^2)),p.nfft,[],[],1/fline);
                figure(p.fig2);
                semilogy(f,abs(pxx)/divisor);
                drawnow
                hold on
                [pxx,f]=nt_spect_plot((x-y)/sqrt(mean(x(:).^2)),p.nfft,[],[],1/fline);
                figure(p.fig2);
                semilogy(f,abs(pxx)/divisor);
                legend('clean', 'removed'); legend boxoff
                set(gca,'ygrid','on','xgrid','on');
                set(gca,'yticklabel',[]); ylabel([]);
                xlabel('frequency (relative to line)');
                yl2=get(gca,'ylim');
                hh=get(gca,'children');
                set(hh(1),'color',[1 .5 .5]); set(hh(2), 'color', [ 0 .7 0]);
                set(hh(2),'linewidth', 2);
                yl(1)=min(yl1(1),yl2(1)); yl(2)=max(yl1(2),yl2(2));
                subplot 121; ylim(yl); subplot 122; ylim(yl);
                drawnow
                
                disp('...done.')
            end
        end
        
        % store in EEG file
        EEG.data = y';
        EEG.etc.zapline.n_remove(i_linefreq,:) = zapline_n_remove_final;
        EEG.etc.zapline.score(i_linefreq,:,:) = scores;
        
        if save_file_on_disk && plot_zapline
            disp('Saving ZapLine figures...')
            
            filenamesplit = strsplit(out_filename,'.set');
            
            savefig(figure(p.fig2),fullfile(out_filepath,[filenamesplit{1}...
                '_' matlab.lang.makeValidName(['zapline_' num2str(linefreq)]) '_spectrum.fig']))
            print(figure(p.fig2),fullfile(out_filepath,[filenamesplit{1}...
                '_' matlab.lang.makeValidName(['zapline_' num2str(linefreq)]) '_spectrum.png']),'-dpng')
            close(p.fig2)
            
            
            savefig(figure(p.fig1),fullfile(out_filepath,[filenamesplit{1}...
                '_' matlab.lang.makeValidName(['zapline_' num2str(linefreq)]) '_scores.fig']))
            print(figure(p.fig1),fullfile(out_filepath,[filenamesplit{1}...
                '_' matlab.lang.makeValidName(['zapline_' num2str(linefreq)]) '_scores.png']),'-dpng')
            close(p.fig1)
            
            disp('...done')
        end
        
    end
    
    % store in EEG file
    EEG.etc.zapline.linefreqs = linefreqs;
    
end

%%

% rename channels if specified
if exist('rename_channels','var') && ~isempty(rename_channels)
    disp('Renaming channels...')
    for i_pair = 1:size(rename_channels,1)
        
        old_chanidx = find(strcmp({EEG.chanlocs.labels},rename_channels{i_pair,1}));
        
        if ~isempty(old_chanidx)
            EEG=pop_chanedit(EEG, 'changefield',{old_chanidx 'labels' rename_channels{i_pair,2}});
        else
            warning(['Did not find channel ' rename_channels{i_pair,1} '. Skipping...'])
        end
        
    end
end

% 1b3) add ref channel as zero if specified
if exist('ref_channel','var') && ~isempty(ref_channel)
    disp('Adding reference channel with zeros...')
    
    EEG.nbchan = EEG.nbchan + 1;
    EEG.data(end+1,:) = zeros(1, EEG.pnts);
    EEG.chanlocs(end+1).labels = ref_channel;
    
    EEG = eeg_checkset(EEG);
    
    disp('...done.')
    
end

% 1c) import chanlocs and copy to urchanlocs
if ~isempty(channel_locations_filepath)
    EEG = pop_chanedit(EEG, 'load',...
        {channel_locations_filepath 'filetype' 'autodetect'});
    disp('Imported channel locations.');
    EEG.urchanlocs = EEG.chanlocs;
else
    eeglab_path = which('eeglab');
    eeglab_path_base = strsplit(eeglab_path,'\eeglab.m');
    standard_channel_locations_path =...
        [eeglab_path_base{1} '\plugins\dipfit2.3\standard_BESA\standard-10-5-cap385.elp'];
    
    EEG=pop_chanedit(EEG,'lookup',standard_channel_locations_path);
end

% this has to happen after loading chanlocs because chanlocs are being completely overwritten in the process
if exist('ref_channel','var')
    
    disp('Declaring ref for all channels...')
    
    [EEG.chanlocs(:).ref] = deal(ref_channel);
    
end


% 1d) change channel types in standard MoBI montage declaring the EOG and EEG channels
for n = 1:length(EEG.chanlocs)
    if ismember(lower(EEG.chanlocs(n).labels), lower(eog_channels))
        EEG.chanlocs(n).type = strcat('EOG');
        disp(['Added channel type: ', EEG.chanlocs(n).labels, ' is EOG electrode!!']);
    elseif exist('ref_channel','var') && strcmpi(EEG.chanlocs(n).labels, ref_channel)
        EEG.chanlocs(n).type = strcat('REF');
        disp(['Added channel type: ', EEG.chanlocs(n).labels, ' is REF electrode!!']);
    else
        EEG.chanlocs(n).type = strcat('EEG');
        disp(['Added channel type: ', EEG.chanlocs(n).labels, ' is EEG electrode.']);
    end
end

EEG = eeg_checkset( EEG );


%% new data set in EEGLAB
[ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, CURRENTSET, 'gui', 'off');
EEG = eeg_checkset( EEG );

% save on disk
if save_file_on_disk
    EEG = pop_saveset( EEG, 'filename',out_filename,'filepath', out_filepath);
    disp('...done')
end

[ALLEEG EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);