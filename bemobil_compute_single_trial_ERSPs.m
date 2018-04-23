% bemobil_compute_single_trial_ERSPs() - Computes single trial ERSP data for a given subject. NO OUTPUT, files will be 
% saved on the disk. 
%
% Usage:
%   >> bemobil_compute_single_trial_ERSPs( input_path , input_filename,  subjects, components_to_use_for_study,...
%       timewarp_name, epochs_info_loadpath, recompute, do_timewarp, dont_warp_but_cut, n_freqs, n_times )
% 
% Inputs:
%   EEG                             - current EEGLAB EEG structure
%   subject                         - subject to be calculated
%   components_to_use_for_study     - which independent components should be used for calculation
%   timewarp_name                   - the name of the used timewarp. this
%                                   will be part of the filepath where the
%                                   ERSPs will be written to
%   epochs_info_loadpath            - filepath to load previously stored epochs information. Useful only if you want to
%                                   store timewarp latencies there, otherwise ignored.
%   recompute                       - boolean, if should recompute, even if
%                                   previous data is present
%   do_timewarp                     - if there should be no timewarping applied at all, set this to false. In this case 
%                                   the timewarp_name is only the output path
%   dont_warp_but_cut               - if a timewarp is present but should only be used to cut the ERSP and leave the
%                                   rest as NAN, set this true
%   n_freqs                         - number of frequencies for ERSP
%   n_times                         - number of timepoints of timewarp
%   fft_options                     - struct with the complete fft options
%                                   (cycles,freqrange,freqscale,padratio,
%                                   alpha,powbase)
%   filepath                        - path where the new single trial ERSP
%                                   base directory will be (output_path =
%                                   [output_filepath 'ERSPs\' timewarp_name
%                                   '\IC_' num2str(IC)];
%
% Outputs:
%
%   NONE
%
% See also:
%   EEGLAB, newtimef, make_timewarp
%
% Authors: Marius Klug, 2018

function bemobil_compute_single_trial_ERSPs( input_path , input_filename,  subjects, components_to_use_for_study, timewarp_name, epochs_info_loadpath, recompute, do_timewarp, dont_warp_but_cut, n_freqs, n_times )



fft_options = struct();
fft_options.cycles = [3 0.5];
fft_options.freqrange = [3 100];
fft_options.freqscale = 'log';
fft_options.padratio = 2;
fft_options.alpha = NaN;
fft_options.powbase = NaN;


if ~exist('ALLEEG','var'); eeglab; end
pop_editoptions( 'option_storedisk', 0, 'option_savetwofiles', 1, 'option_saveversion6', 0, 'option_single', 0, 'option_memmapdata', 0, 'option_eegobject', 0, 'option_computeica', 1, 'option_scaleicarms', 1, 'option_rememberfolder', 1, 'option_donotusetoolboxes', 0, 'option_checkversion', 1, 'option_chat', 1);

if do_timewarp
    try
        load([input_path timewarp_name],'timeWarp')
    catch
       error(['timeWarp struct could not be loaded using ''' input_path timewarp_name '''!']) 
    end
    try
        load([input_path timewarp_name '_latencyMeans'],'latencyMeans')
    catch
        error(['latencyMeans struct could not be loaded using ''' input_path timewarp_name '''!'])
    end
else
end

try 
    % load epoch_info. load stores into a struct, so the first element of the struct has to be taken
    epochs_info = load(epochs_info_loadpath);
    epoch_info_fields = fieldnames(epochs_info);
    epochs_info = epochs_info.(epoch_info_fields{1});
    
    epochs_info_present = true;
catch
    warning('Loading epoch info failed, no timewarp info will be saved in epochs!')
    epochs_info_present = false;
end

for subject = subjects
    STUDY = []; CURRENTSTUDY = 0; ALLEEG = []; EEG=[]; CURRENTSET=[];
    
    filepath = [input_path num2str(subject) '\'];
    
    EEG = pop_loadset('filename',input_filename,'filepath',filepath);
    
    % compute newtimef data
    for IC = components_to_use_for_study
        tic % start checking the time
        disp(['Subject: ' num2str(subject)])
        disp(['IC: ' num2str(IC)])
        output_path = [filepath 'ERSPs\' timewarp_name '\IC_' num2str(IC)];
        directory_content = dir(output_path);
        
        % check if this is already calculated
        already_calculated = false;
        for data = 1:length(directory_content)
            if strcmp(directory_content(data).name,'all_epochs_ersp.mat')
                already_calculated = true;
                break
            end
        end
        
        if already_calculated && ~recompute
            continue
        end
        
        all_epochs_ersp = nan(length(EEG.epoch),n_freqs,n_times);
        
        if do_timewarp
            this_subject_timewarp_latencies = timeWarp(subject).latencies;

            % for some reason this is necessary, maybe a bug in maketimewarp, anyways, there exist instances, where the
            % last timewarp marker is 1 frame after the end of the epoch...
            if this_subject_timewarp_latencies(this_subject_timewarp_latencies < EEG.times(1))
                warning('Some timewarp markers were before epoch start.  Corrected to epoch limits...')
                disp('Were:')
                this_subject_timewarp_latencies(this_subject_timewarp_latencies < EEG.times(1))
                disp('Are:')
                EEG.times(1)+(1/EEG.srate*1000*2)
                this_subject_timewarp_latencies(this_subject_timewarp_latencies < EEG.times(1)) = EEG.times(1)+(1/EEG.srate*1000*2);
            end
            if this_subject_timewarp_latencies(this_subject_timewarp_latencies > EEG.times(end))
                warning('Some timewarp markers were after epoch end. Corrected to epoch limits...')
                disp('Were:')
                this_subject_timewarp_latencies(this_subject_timewarp_latencies > EEG.times(end))
                disp('Are:')
                EEG.times(end)-(1/EEG.srate*1000*2)
                this_subject_timewarp_latencies(this_subject_timewarp_latencies > EEG.times(end)) = EEG.times(end)-(1/EEG.srate*1000*2);
            end

            % not all epochs do have a timewarp, since there might have been a wrong marker order or nonexistent markers
            for epoch_index = 1:length(timeWarp(subject).epochs) 

                fprintf('.')

                epoch = timeWarp(subject).epochs(epoch_index); % find the actual epoch number

                thisTimeWarp = this_subject_timewarp_latencies(epoch_index,:); % this is the timewarp latency for this epoch

                if dont_warp_but_cut
                    
                    % do the timefreq analysis without timewarp
                    [full_ersp,~,~,times,~,~,~] = newtimef(EEG.icaact(IC,:,epoch),...
                        EEG.pnts,...
                        [EEG.times(1) EEG.times(end)],...
                        EEG.srate,...
                        'cycles',fft_options.cycles,...
                        'freqs',fft_options.freqrange,...
                        'nfreqs',n_freqs,...
                        'timesout',n_times,...
                        'freqscale',fft_options.freqscale,...
                        'padratio',fft_options.padratio,...
                        'alpha',fft_options.alpha,...
                        'powbase',fft_options.powbase,...
                        'baseline',NaN,... % no baseline, since that is only a subtraction of the freq values, we do it manually
                        'plotersp','off',...
                        'plotitc','off',...
                        'verbose','off');
                    
                    % have a NaN ERSP
                    ersp = NaN(n_freqs,n_times);
                    
                    % fill it with the ersp data up to the last timewarp latency frame
                    ersp(:,1:find(times>=thisTimeWarp(end),1,'first')) = full_ersp(:,1:find(times>=thisTimeWarp(end),1,'first'));
                    
                else
                    
                    % do normal timewarp
                    [ersp,~,~,~,~,~,~] = newtimef(EEG.icaact(IC,:,epoch),...
                        EEG.pnts,...
                        [EEG.times(1) EEG.times(end)],...
                        EEG.srate,...
                        'cycles',fft_options.cycles,...
                        'freqs',fft_options.freqrange,...
                        'nfreqs',n_freqs,...
                        'timewarp',thisTimeWarp,...
                        'timesout',n_times,...
                        'freqscale',fft_options.freqscale,...
                        'padratio',fft_options.padratio,...
                        'alpha',fft_options.alpha,...
                        'powbase',fft_options.powbase,...
                        'timewarpms',latencyMeans,... % mean latencies of the timewarp over all epochs (calculated above)
                        'baseline',NaN,... % no baseline, since that is only a subtraction of the freq values, we do it manually
                        'plotersp','off',...
                        'plotitc','off',...
                        'verbose','off');
                end

                % store ersp in the data set that will be saved
                all_epochs_ersp(timeWarp(subject).epochs(epoch_index),:,:) = single(ersp); % single saves disk space
                
                if epochs_info_present
                    % store timewarp latency infos
                    epochs_info(timeWarp(subject).epochs(epoch_index)).timeWarpLatencies = thisTimeWarp;
                end

            end
        else
            
            for epoch = 1:length(EEG.epoch)

                fprintf('.')

                % do timefreq analysis without timewarp
                [ersp,~,~,~,~,~,~] = newtimef(EEG.icaact(IC,:,epoch),...
                    EEG.pnts,...
                    [EEG.times(1) EEG.times(end)],...
                    EEG.srate,...
                    'cycles',fft_options.cycles,...
                    'freqs',fft_options.freqrange,...
                    'nfreqs',n_freqs,...
                    'timesout',n_times,...
                    'freqscale',fft_options.freqscale,...
                    'padratio',fft_options.padratio,...
                    'alpha',fft_options.alpha,...
                    'powbase',fft_options.powbase,...
                    'baseline',NaN,... % no baseline, since that is only a subtraction of the freq values, we do it manually
                    'plotersp','off',...
                    'plotitc','off',...
                    'verbose','off');

                % store ersp in the data set that will be saved
                all_epochs_ersp(epoch,:,:) = single(ersp); % single saves disk space
                
                if epochs_info_present
                    
                    % no timewarp was done, store that
                    epochs_info(epoch).timeWarpLatencies = [];
                end

            end
        end
        
        mkdir(output_path)
        save([output_path '\all_epochs_ersp'], 'all_epochs_ersp');
        
        lastduration = toc; % stop checking the time and plot estimated time of arrival
        eta = lastduration * ( length(subjects)*length(components_to_use_for_study)- (((find(subject==subjects)-1)*length(components_to_use_for_study))+IC));
        disp('\nLast duration:')
        disp(lastduration)
        disp('ETA (h):')
        disp(eta/3600)
    end
end

% those are always the same and only need to be saved once
save([input_path '\times'], 'times');
save([input_path '\freqs'], 'freqs');
disp('Done.')