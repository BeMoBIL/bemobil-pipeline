% bemobil_compute_single_trial_ERSPs_with_timewarp() - Computes single
% trial ERSP data for a given subject. NO OUTPUT, files will be saved on
% the disk.
%
% Usage:
%   >> bemobil_compute_single_trial_ERSPs_with_timewarp( EEG , subject, components_to_use_for_study, timeWarp, timewarp_name, latencyMeans, recompute, n_freqs, n_times, fft_options, filepath )
%
% Inputs:
%   EEG                             - current EEGLAB EEG structure
%   subject                         - subject to be calculated
%   components_to_use_for_study     - which independent components should be used for calculation
%   timewarp_name                   - the name of the used timewarp. this
%                                   will be part of the filepath where the
%                                   ERSPs will be written to
%   recompute                       - boolean, if should recompute, even if
%                                   previous data is present
%   n_freqs                         - number of frequencies for ERSP
%   n_times                         - number of timepoints of timewarp
%   fft_options                     - struct with the complete fft options
%                                   (cycles,freqrange,freqscale,padratio,
%                                   alpha,powbase)
%   output_filepath                 - path where the new single trial ERSP
%                                   base directory will be (output_path =
%                                   [output_filepath 'ERSPs\' timewarp_name
%                                   '\IC_' num2str(IC)];
%
% Outputs:
%
%   NONE
%
% See also:
%   EEGLAB, newtimef
%
% Authors: Marius Klug, 2018

function bemobil_compute_single_trial_ERSPs_with_timewarp( input_path , input_filename,  subjects, components_to_use_for_study, timewarp_name, recompute, n_freqs, n_times )

fft_options = struct();
fft_options.cycles = [3 0.5];
fft_options.freqrange = [3 100];
fft_options.freqscale = 'log';
fft_options.padratio = 2;
fft_options.alpha = NaN;
fft_options.powbase = NaN;

eeglab

load([input_path timewarp_name],'timeWarp')
load([input_path timewarp_name '_latencyMeans'],'latencyMeans')

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
        
        for epoch_index = 1:length(timeWarp(subject).epochs) % not all epochs do have a timewarp, since there might have been a wrong marker order or nonexistent markers
            
            fprintf('.')
            
            epoch = timeWarp(subject).epochs(epoch_index); % find the actual epoch number
            
            thisTimeWarp = timeWarp(subject).latencies(epoch_index,:); % this is the timewarp latency for this epoch
            
            %             output_path = [filepath 'ERSPs\' timewarp_name_1 '\IC_' num2str(IC) '\epoch_' num2str(epoch)];
            %
            %             mkdir(output_path);
            %
            %             if size(dir(output_path),1) == 4
            %                 % 4 is correct for current script, so this epoch can be skipped
            %                 continue;
            %             end
            %
            
            [ersp,itc,powbase,times,freqs,erspboot,itcboot] = newtimef(EEG.icaact(IC,:,epoch),...
                EEG.pnts,...
                [-1000 13000],...
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
            
            all_epochs_ersp(timeWarp(subject).epochs(epoch_index),:,:) = single(ersp); % single saves disk space
            
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