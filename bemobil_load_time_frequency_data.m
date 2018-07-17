function time_frequency_data = bemobil_load_time_frequency_data(input_path_epochs, nepochs, freqs, times, n_freqs, n_times, timewarp_name, trial_normalization, baseline_start_end, experiment_conditions_to_test, epoch_rejections_this_subject, epoch_rejections_this_subject_for_baseline, subject, IC)

% load ERSPs and compute study ERSPs
% load('P:\Lukas_Gehrke\studies\Spot_Rotation\data\4_single_subject_analysis\ERSPs\ERSPbaselines','ERSPbaselines');
load([input_path_epochs timewarp_name '_latencyMeans'],'latencyMeans');

load([input_path_epochs timewarp_name],'timeWarp')

load([input_path_epochs timewarp_name '_timeWarpOutlierEpoch'],'timeWarpOutlierEpoch')
% load([input_path_epochs timewarp_name '_outlierEpoch'],'outlierEpoch')
load([input_path_epochs timewarp_name '_subjectLatencyMeans'],'subjectLatencyMeans')
load([input_path_epochs timewarp_name '_latencyMeans'],'latencyMeans')

time_frequency_data = struct('erspboot',[]);

% clear cluster_ersp_across_subjects_power
% clear subject_ersps_across_ICs_power
% clear baseline
% clear residual_variances

% baseline_to_use = 'grand average';
% baseline_to_use = 'trial';
% baseline_to_use = 'resting condition';
% baseline_to_use = 'none';
% baseline_to_use = 'trial normalization grand average prestim';

% experiment_conditions_to_test = {'body','joystick'};

% residual_variance_threshold = 0.2;


% for cluster = clusters
%
%     disp(['Cluster ' num2str(cluster) ' of [' num2str(clusters) ']'])
%
%     ERSP_subjects = clustering_solution(cluster).sets(ismember(clustering_solution(cluster).sets, subjects));
%     ICs = clustering_solution(cluster).comps(ismember(clustering_solution(cluster).sets, subjects));
%     old_subject = 0;
%
%     subject_ersps_across_epochs_power = NaN(nsubjects,maxIC,length(experiment_conditions_to_test),n_freqs,n_times);
%     subject_ersp_thisIC_all_epochs_power = NaN(nepochs,length(experiment_conditions_to_test),n_freqs,n_times);
%
%     ICs_thrown_out_this_cluster = 0;
%
%     indIC = 1;
%     while indIC <= length(ICs)
%
%         % current IC
%         IC = ICs(indIC);
%
%         % set of the current IC
%         subject = ERSP_subjects(indIC);
%
%         % RV of current IC
%         if ~isempty(ALLEEG)
%
%             residual_variance = ALLEEG(subject).dipfit.model(IC).rv;
%
%         else
%
%             disp('RV not taken into account. Load the whole study to enable this.')
%
%             residual_variance = 0;
%
%
%         end
%         residual_variances(cluster,indIC) = residual_variance;
%
%         disp(['IC #' num2str(indIC) ' of ' num2str(length(ICs)) ', RV: ' num2str(round(residual_variance*100)) '%'])
%
%         if residual_variance <= residual_variance_threshold

filepath = [input_path_epochs num2str(subject) '\'];

subject_ersp_thisIC_all_epochs_power = NaN(nepochs,max(1,length(experiment_conditions_to_test)),n_freqs,n_times);
%
%             % THIS IS THE PRECALCULATED BASELINE FOR THE BODY(1)/JOYSTICK(2) CONDITION!!!
%             for condition = 1:size(ERSPbaselines,2)
%
%                 baseline(:,condition) = squeeze(ERSPbaselines(subject,condition,:));
%                 % baselines are in power values already
%
%             end
%
%             baseline = 10.^(baseline/10);

for epoch_index = 1:nepochs
    
    if ~timeWarpOutlierEpoch(subject,epoch_index)
        
        input_path = [filepath 'ERSPs\' timewarp_name '\IC_' num2str(IC) '\epoch_' num2str(epoch_index)];
        
        % load epoch info, but not all epochs are present, due to timewarping mismatches
        try
            load([input_path '\epoch_info'], 'epoch_info')
            load([input_path '\ersp'], 'ersp');
            
            
            % transform from dB to power values
            power = 10.^(ersp/10);
            
            if trial_normalization
                full_trial_baseline = mean(power,2);
                power = power ./ full_trial_baseline;
            end
            
            if ~isempty(experiment_conditions_to_test)
                for condition = 1:length(experiment_conditions_to_test)
                    
                    if strcmp(epoch_info.eventexperiment_condition,experiment_conditions_to_test{condition})
                        subject_ersp_thisIC_all_epochs_power(epoch_index,condition,:,:)=power;
                    end
                    
                end
            else
                subject_ersp_thisIC_all_epochs_power(epoch_index,1,:,:)=power;
            end
            
            
        catch
            warning('There was an error in loading or processing an epoch!')
        end
    end
    
end

% this is the mean ersp of this IC by this subject across all epochs that satisfy the conditions and
% are suitable for the ERSP
subject_ersps_across_epochs_power = squeeze(nanmean(subject_ersp_thisIC_all_epochs_power(~logical(epoch_rejections_this_subject),:,:,:),1));

% this is the mean ersp of this IC by this subject across all epochs that satisfy the conditions and
% are suitable for the baseline of the ERSP
subject_ersps_across_epochs_power_for_baseline = squeeze(nanmean(subject_ersp_thisIC_all_epochs_power(~logical(epoch_rejections_this_subject_for_baseline),:,:,:),1));

if isempty(experiment_conditions_to_test)
    time_frequency_data.powbase = squeeze(mean(subject_ersps_across_epochs_power_for_baseline(:,find(times>baseline_start_end(1),1,'first'):find(times<=baseline_start_end(2),1,'last')),2));
    subject_baseline_corrected_ersps_across_epochs_power = subject_ersps_across_epochs_power ./ time_frequency_data.powbase;
    time_frequency_data.ersp = single(10.*log10(subject_baseline_corrected_ersps_across_epochs_power));
    
else
    baselines_power = squeeze(mean(subject_ersps_across_epochs_power(:,:,find(times>baseline_start_end(1),1,'first'):find(times<=baseline_start_end(2),1,'last')),3));
end

time_frequency_data.freqs = freqs;
time_frequency_data.times = times;


%         else
%
%             disp(['IC not taken into account due to RV threshold of ' num2str(residual_variance_threshold*100) '%.'])
%             ICs_thrown_out_this_cluster = ICs_thrown_out_this_cluster + 1;
%
%         end
%
%         % increase counter
%         indIC = indIC+1;
%     end
%
%     subject_ersps_across_ICs_power(cluster,:,:,:,:) = squeeze(nanmean(subject_ersps_across_epochs_power,2));
%
%     %     if length(subjects) > 1 || subjects ~= 1
%
%     cluster_ersp_across_subjects_power = squeeze(nanmean(subject_ersps_across_ICs_power,2));
%
%     disp(['Average RV in this cluster without thresholding: ' num2str(round(nanmean(squeeze(residual_variances(cluster,:)))*100)) '%']);
%     disp(['ICs thrown out: ' num2str(ICs_thrown_out_this_cluster) '/' num2str(length(ICs)) ' (' num2str(round(ICs_thrown_out_this_cluster/length(ICs)*100)) '%)'])
%
%
%     else
%
%         data_to_plot_power = subject_ersps_across_ICs_power;
%
%         if strcmp(baseline_to_use,'trial')
%
%             data_to_plot =  10.*log10(data_to_plot_power);
%
%         elseif strcmp(baseline_to_use,'grand average')
%
%             baseline_corrected = data_to_plot_power ./ mean(data_to_plot_power(:,1:find(times<0,1,'last')),2);
%             data_to_plot =  10.*log10(baseline_corrected);
%
%         end
%
%     end;



%     startIndex = 1;
%     endIndex = find(times<latencyMeans(3),1,'last');
%     scale_max = max(max(abs(data_to_plot(:,startIndex:endIndex))))/2;

%     figure;
%     imagesclogy(times(startIndex:endIndex),freqs,data_to_plot(:,startIndex:endIndex),[-scale_max scale_max]);
%     axis xy;
%     vline(latencyMeans,'black');title(['Cluster: ' num2str(cluster) ', baseline: ' baseline_to_use ', condition: ' condition_to_plot ', symmetric scaling, ERSP subjects: ' num2str(unique(ERSP_subjects))]);
%     cbar;

%     figure;
%     imagesclogy(times(startIndex:endIndex),freqs,data_to_plot(:,startIndex:endIndex));
%     axis xy;
%     vline(latencyMeans,'black');title(['Cluster: ' num2str(cluster) ', baseline: ' baseline_to_use ', condition: ' condition_to_plot ', auto scaling, ERSP subjects: ' num2str(unique(ERSP_subjects))]);
%     cbar;

% end


end