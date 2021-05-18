% bemobil_compute_single_trial_ERSPs() - Computes single trial ERSP data for a given subject. NO OUTPUT, files will be 
% saved on the disk. 
%
% Usage:
% function bemobil_compute_single_trial_ERSPs( input_path , input_filename,  subjects, components_to_use_for_study,...
%     output_foldername, timewarp_latency_loadpath, epochs_info_filename_input, epochs_info_filename_output, recompute, do_timewarp,...
%     dont_warp_but_cut, n_freqs, n_times )
% 
% Inputs:
%   input_path                      - path to the EEG data sets (without the subject number)
%   input_filename                  - name of the epoch data set
%   subjects                        - vector of subjects that should be calculated
%   components_to_use_for_study     - which independent components should be used for calculation
%   output_foldername               - name of the folder where the ERSPs will be saved 
%   timewarp_latency_loadpath       - the full path to the timewarp latencies INCLUDING the name
%   epochs_info_filename_input      - filename of previously stored epochs information. Useful only if you want to
%                                   store timewarp latencies there, otherwise ignored. 
%   epochs_info_filename_output     - filename to save new epochs information with timewarp info per epoch. Useful only 
%                                   if you want to store timewarp latencies there, otherwise ignored. 
%   recompute                       - boolean, if should recompute, even if
%                                   previous data is present
%   has_timewarp_latencies          - set true if there are precomputed timewarp latencies present. Will lead to
%                                   timewarped ERSPs
%   dont_warp_but_cut               - if a timewarp is present but should only be used to cut the ERSP and leave the
%                                   rest as NAN, set this true
%   n_freqs                         - number of frequencies for ERSP
%   n_times                         - number of timepoints of timewarp
%
% Outputs:
%
%   NONE!
%   files are STORED on the DISK.
%   output_path = [input_path num2str(subject) '\ERSPs\' output_foldername '\IC_' num2str(IC)];
%   save([input_path '\' num2str(subject) '\' epochs_info_filename_output], 'epochs_info');
%
% See also:
%   EEGLAB, newtimef, make_timewarp
%
% Authors: Marius Klug, 2018

function bemobil_compute_single_trial_ERSPs( input_path , input_filename,  subjects, components_to_use_for_study,...
    output_foldername, timewarp_latency_loadpath, epochs_info_filename_input, epochs_info_filename_output, recompute, has_timewarp_latencies,...
    dont_warp_but_cut, n_freqs, n_times )



fft_options = struct();
fft_options.cycles = [3 0.5];
fft_options.freqrange = [3 100];
fft_options.freqscale = 'log';
fft_options.padratio = 2;
fft_options.alpha = NaN;
fft_options.powbase = NaN;


if ~exist('ALLEEG','var'); eeglab; end
pop_editoptions( 'option_storedisk', 0, 'option_savetwofiles', 1, 'option_saveversion6', 0, 'option_single', 0, 'option_memmapdata', 0, 'option_eegobject', 0, 'option_computeica', 1, 'option_scaleicarms', 1, 'option_rememberfolder', 1, 'option_donotusetoolboxes', 0, 'option_checkversion', 1, 'option_chat', 1);

if has_timewarp_latencies
    try
        load(timewarp_latency_loadpath,'timeWarp')
    catch
       error(['timeWarp struct could not be loaded using ''' timewarp_latency_loadpath '''!']) 
    end
    try
        load([timewarp_latency_loadpath '_latencyMeans'],'latencyMeans')
    catch
        error(['latencyMeans struct could not be loaded using ''' timewarp_latency_loadpath '''!'])
    end
else
end

for subject = subjects
    STUDY = []; CURRENTSTUDY = 0; ALLEEG = []; EEG=[]; CURRENTSET=[];
    
    filepath = [input_path num2str(subject) '\'];
    
    EEG = pop_loadset('filename',input_filename,'filepath',filepath);
    
    try 
        % load epoch_info. load stores into a struct, so the first element of the struct has to be taken
        epochs_info = load([input_path '\' num2str(subject) '\' epochs_info_filename_input]);
        epoch_info_fields = fieldnames(epochs_info);
        epochs_info = epochs_info.(epoch_info_fields{1});

        epochs_info_present = true;
    catch
        warning('Loading epoch info failed, no timewarp info will be saved in epochs!')
        epochs_info_present = false;
    end
    
    % compute newtimef data
    for IC = components_to_use_for_study
        tic % start checking the time
        disp(['Subject: ' num2str(subject)])
        disp(['IC: ' num2str(IC)])
        
        output_path = [filepath 'ERSPs\' output_foldername '\IC_' num2str(IC)];
        
        % check if this is already calculated
        directory_content = dir(output_path);
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
        
        if has_timewarp_latencies
            this_subject_timewarp_latencies = timeWarp(subject).latencies;

            % for some reason this is necessary, maybe a bug in maketimewarp, anyways, there exist instances, where the
            % last timewarp marker is 1 frame after the end of the epoch...
            
            minTime = EEG.times(1) + (2 * (1 / EEG.srate) * 1000);
            maxTime = EEG.times(end) - (2 * (1 / EEG.srate) * 1000);
            
            if this_subject_timewarp_latencies(this_subject_timewarp_latencies < minTime)
                warning('Some timewarp markers were before epoch start.  Corrected to epoch limits...')
                disp('Were:')
                this_subject_timewarp_latencies(this_subject_timewarp_latencies < minTime)
                disp(['Are: ' num2str(minTime)])
                this_subject_timewarp_latencies(this_subject_timewarp_latencies < minTime) = minTime;
            end
            if this_subject_timewarp_latencies(this_subject_timewarp_latencies > maxTime)
                warning('Some timewarp markers were after epoch end. Corrected to epoch limits...')
                disp('Were:')
                this_subject_timewarp_latencies(this_subject_timewarp_latencies > maxTime)
                disp(['Are: ' num2str(maxTime)])
                this_subject_timewarp_latencies(this_subject_timewarp_latencies > maxTime) = maxTime;
            end

            % not all epochs do have a timewarp, since there might have been a wrong marker order or nonexistent markers
            for epoch_index = 1:length(timeWarp(subject).epochs) 

                fprintf('.')

                epoch = timeWarp(subject).epochs(epoch_index); % find the actual epoch number

                thisTimeWarp = this_subject_timewarp_latencies(epoch_index,:); % this is the timewarp latency for this epoch

                if dont_warp_but_cut
                    
                    % do the timefreq analysis without timewarp
                    [full_ersp,~,~,times,freqs,~,~] = newtimef(EEG.icaact(IC,:,epoch),...
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
                    
                    assert(all(size(thisTimeWarp)==size(latencyMeans)),'The timewarp matrix of this epoch has a different size than the latency means to warp.')
                    
                    % do normal timewarp
                    [ersp,~,~,times,freqs,~,~] = newtimef(EEG.icaact(IC,:,epoch),...
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
                [ersp,~,~,times,freqs,~,~] = newtimef(EEG.icaact(IC,:,epoch),...
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
    
    % this is a little dirty, since the epoch info will be filled each IC again, but it is identical, so won't matter
    if epochs_info_present
        save([input_path '\' num2str(subject) '\' epochs_info_filename_output], 'epochs_info');
    end
    
end

% those are always the same and only need to be saved once
try
    save([input_path '\times'], 'times');
    save([input_path '\freqs'], 'freqs');
catch
    warning('No trials have been actually computed, nothing was saved on the disk. Probably everything had been computed before and recompute was set to 0.')
end
disp('Done.')