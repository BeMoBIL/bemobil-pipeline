function time_frequency_data = bemobil_load_time_frequency_data_single_subject(input_path_epochs, nepochs, freqs, times, timewarp_name, trial_normalization, baseline_start_end, experiment_conditions_to_test, epoch_rejections_this_subject, epoch_rejections_this_subject_for_baseline, subject, IC)

% load ERSPs and compute study ERSPs
% load('P:\Lukas_Gehrke\studies\Spot_Rotation\data\4_single_subject_analysis\ERSPs\ERSPbaselines','ERSPbaselines');
load([input_path_epochs timewarp_name '_latencyMeans'],'latencyMeans');

load([input_path_epochs timewarp_name],'timeWarp')

load([input_path_epochs timewarp_name '_timeWarpOutlierEpoch'],'timeWarpOutlierEpoch')
load([input_path_epochs timewarp_name '_subjectLatencyMeans'],'subjectLatencyMeans')
load([input_path_epochs timewarp_name '_latencyMeans'],'latencyMeans')

time_frequency_data.baseline_start_end = baseline_start_end;
time_frequency_data.epoch_rejections_for_ERSPs = epoch_rejections_this_subject;
time_frequency_data.epoch_rejections_for_baseline = epoch_rejections_this_subject_for_baseline;
time_frequency_data.timewarp_name = timewarp_name;

time_frequency_data.freqs = freqs;
time_frequency_data.times = times;
time_frequency_data.experiment_conditions_to_test = experiment_conditions_to_test;
time_frequency_data.grand_average.erspboot = [];

filepath = [input_path_epochs num2str(subject) '\'];

subject_ersp_thisIC_all_epochs_power = NaN(nepochs,length(freqs),length(times));

experiment_conditions = zeros(1,280);

% baseline_to_use = 'grand average';
% baseline_to_use = 'trial';
% baseline_to_use = 'resting condition';
% baseline_to_use = 'none';
% baseline_to_use = 'trial normalization grand average prestim';

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
            
            subject_ersp_thisIC_all_epochs_power(epoch_index,:,:)=power;
            
            if ~isempty(experiment_conditions_to_test)
                for condition = 1:length(experiment_conditions_to_test)
                    if strcmp(epoch_info.eventexperiment_condition,experiment_conditions_to_test{condition})
                        experiment_conditions(epoch_index) = condition;
                    end
                    
                end
            else
                
            end
            
            
        catch
            warning('There was an error in loading or processing an epoch!')
        end
    end
    
end

% this is the mean ersp of this IC by this subject across all epochs that
% are suitable for the ERSP
subject_ersps_across_epochs_power_grand_average = squeeze(nanmean(subject_ersp_thisIC_all_epochs_power(~logical(epoch_rejections_this_subject),:,:),1));

% this is the mean ersp of this IC by this subject across all epochs that
% are suitable for the baseline of the ERSP
subject_ersps_across_epochs_power_grand_average_for_baseline = squeeze(nanmean(subject_ersp_thisIC_all_epochs_power(~logical(epoch_rejections_this_subject_for_baseline),:,:),1));


time_frequency_data.grand_average.base_power = single(squeeze(mean(subject_ersps_across_epochs_power_grand_average_for_baseline(:,find(times>baseline_start_end(1),1,'first'):find(times<=baseline_start_end(2),1,'last')),2)));
time_frequency_data.grand_average.base_power_ersp = single(repmat(time_frequency_data.grand_average.base_power,1,length(times)));
time_frequency_data.grand_average.base_dB_ersp = single(10.*log10(time_frequency_data.grand_average.base_power_ersp)); % this is necessary for bootstrapping
time_frequency_data.grand_average.raw_power = single(subject_ersps_across_epochs_power_grand_average);
time_frequency_data.grand_average.raw_dB = single(10.*log10(subject_ersps_across_epochs_power_grand_average));
time_frequency_data.grand_average.ersp_power = single(subject_ersps_across_epochs_power_grand_average ./ time_frequency_data.grand_average.base_power);
time_frequency_data.grand_average.ersp = single(10.*log10(time_frequency_data.grand_average.ersp_power));
time_frequency_data.grand_average.condition = 'grand average';


for condition = 1:length(experiment_conditions_to_test)
    
    % this is the mean ersp of this IC by this subject across all epochs that
    % are suitable for the ERSP
    subject_ersps_across_epochs_power_this_condition = squeeze(nanmean(subject_ersp_thisIC_all_epochs_power(~logical(epoch_rejections_this_subject)&experiment_conditions==condition,:,:,:),1));
    
    % this is the mean ersp of this IC by this subject across all epochs that
    % are suitable for the baseline of the ERSP
    subject_ersps_across_epochs_power_this_for_baseline = squeeze(nanmean(subject_ersp_thisIC_all_epochs_power(~logical(epoch_rejections_this_subject_for_baseline)&experiment_conditions==condition,:,:,:),1));
    
    time_frequency_data.(['condition_' num2str(condition)]).base_power = single(squeeze(mean(subject_ersps_across_epochs_power_this_for_baseline(:,find(times>baseline_start_end(1),1,'first'):find(times<=baseline_start_end(2),1,'last')),2)));
    time_frequency_data.(['condition_' num2str(condition)]).base_power_ersp = single(repmat(time_frequency_data.(['condition_' num2str(condition)]).base_power,1,length(times)));
    time_frequency_data.(['condition_' num2str(condition)]).base_dB_ersp = single(10.*log10(time_frequency_data.(['condition_' num2str(condition)]).base_power_ersp)); % this is necessary for bootstrapping
    time_frequency_data.(['condition_' num2str(condition)]).raw_power = single(subject_ersps_across_epochs_power_this_condition);
    time_frequency_data.(['condition_' num2str(condition)]).raw_dB = single(10.*log10(subject_ersps_across_epochs_power_this_condition));
    time_frequency_data.(['condition_' num2str(condition)]).ersp_power = single(subject_ersps_across_epochs_power_this_condition ./ time_frequency_data.(['condition_' num2str(condition)]).base_power);
    time_frequency_data.(['condition_' num2str(condition)]).ersp_with_grand_average_baseline_power = single(subject_ersps_across_epochs_power_this_condition ./ time_frequency_data.grand_average.base_power);
    time_frequency_data.(['condition_' num2str(condition)]).ersp = single(10.*log10(time_frequency_data.(['condition_' num2str(condition)]).ersp_power));
    time_frequency_data.(['condition_' num2str(condition)]).ersp_with_grand_average_baseline = single(10.*log10(time_frequency_data.(['condition_' num2str(condition)]).ersp_with_grand_average_baseline_power));
    time_frequency_data.(['condition_' num2str(condition)]).condition = experiment_conditions_to_test{condition};
end

% difference plots
combinations = nchoosek(experiment_conditions_to_test,2);
combinations_indices = nchoosek(1:length(experiment_conditions_to_test),2);
time_frequency_data.differences.combinations = combinations;
for combination = 1:size(combinations,1)
%     time_frequency_data.differences.(['combination_' num2str(combination)]).ersp = single(time_frequency_data.(['condition_' num2str(combinations_indices(combination,1))]).ersp - time_frequency_data.(['condition_' num2str(combinations_indices(combination,2))]).ersp);
%     time_frequency_data.differences.(['combination_' num2str(combination)]).ersp_power = single(10.^(time_frequency_data.differences.(['combination_' num2str(combination)]).ersp/10)); % power = 10.^(dB/10)
%     time_frequency_data.differences.(['combination_' num2str(combination)]).base_power = single((time_frequency_data.(['condition_' num2str(combinations_indices(combination,1))]).base_power + time_frequency_data.(['condition_' num2str(combinations_indices(combination,2))]).base_power)/2);
%     time_frequency_data.differences.(['combination_' num2str(combination)]).base_power_ersp = single(repmat(time_frequency_data.differences.(['combination_' num2str(combination)]).base_power,1,length(times)));
%     time_frequency_data.differences.(['combination_' num2str(combination)]).base_dB_ersp = single(10.*log10(time_frequency_data.differences.(['combination_' num2str(combination)]).base_power_ersp)); % this is necessary for bootstrapping
%     time_frequency_data.differences.(['combination_' num2str(combination)]).raw_power = single(time_frequency_data.differences.(['combination_' num2str(combination)]).ersp_power .* time_frequency_data.differences.(['combination_' num2str(combination)]).base_power);
%     time_frequency_data.differences.(['combination_' num2str(combination)]).raw_dB = single(10.*log10(time_frequency_data.differences.(['combination_' num2str(combination)]).raw_power)); % dB = 10.*log10(power)
%     time_frequency_data.differences.(['combination_' num2str(combination)]).condition = [combinations{combination,1} ' - ' combinations{combination,2}];
    
    time_frequency_data.differences.(['combination_' num2str(combination)]).ersp_power = single(time_frequency_data.(['condition_' num2str(combinations_indices(combination,1))]).ersp_power ./ time_frequency_data.(['condition_' num2str(combinations_indices(combination,2))]).ersp_power);
    
    time_frequency_data.differences.(['combination_' num2str(combination)]).ersp = single(10.*log10(time_frequency_data.differences.(['combination_' num2str(combination)]).ersp_power));
    
    time_frequency_data.differences.(['combination_' num2str(combination)]).base_power = single(squeeze(mean(time_frequency_data.differences.(['combination_' num2str(combination)]).ersp_power(:,find(times>baseline_start_end(1),1,'first'):find(times<=baseline_start_end(2),1,'last')),2)));
    time_frequency_data.differences.(['combination_' num2str(combination)]).base_power_ersp = single(repmat(time_frequency_data.differences.(['combination_' num2str(combination)]).base_power,1,length(times)));
    time_frequency_data.differences.(['combination_' num2str(combination)]).base_dB_ersp = single(10.*log10(time_frequency_data.differences.(['combination_' num2str(combination)]).base_power_ersp)); % this is necessary for bootstrapping
    time_frequency_data.differences.(['combination_' num2str(combination)]).raw_power = single(time_frequency_data.differences.(['combination_' num2str(combination)]).ersp_power .* time_frequency_data.differences.(['combination_' num2str(combination)]).base_power);
    time_frequency_data.differences.(['combination_' num2str(combination)]).raw_dB = single(10.*log10(time_frequency_data.differences.(['combination_' num2str(combination)]).raw_power)); % dB = 10.*log10(power)
    time_frequency_data.differences.(['combination_' num2str(combination)]).condition = [combinations{combination,1} ' - ' combinations{combination,2}];
end

% statistics


% % create full ERSP plots of pure baselines and take the uncorrected ERSP
% all_subjects_prestimulus_baseline_ERSP_dB = repmat(squeeze(all_subjects_prestimulus_baseline_dB(cluster,:,:,:)),[1, 1, 1, length(startTimeIndex:endTimeIndex)]);
% all_subjects_ERSP_uncorrected_ERSP_dB_this_cluster = squeeze(all_subjects_ERSP_uncorrected_dB(cluster,:,:,:,:));
% all_subjects_ERSP_baseline_corrected_dB_this_cluster = squeeze(all_subjects_ERSP_baseline_corrected_dB(cluster,:,:,:,:));
% 
% % permute data such that the subject is the last dimension (necessary for statcond() )
% all_subjects_prestimulus_baseline_ERSP_dB = permute(all_subjects_prestimulus_baseline_ERSP_dB,[1 3 4 2]);
% all_subjects_ERSP_uncorrected_ERSP_dB_this_cluster = permute(all_subjects_ERSP_uncorrected_ERSP_dB_this_cluster,[1 3 4 2]);
% all_subjects_ERSP_baseline_corrected_dB_this_cluster = permute(all_subjects_ERSP_baseline_corrected_dB_this_cluster,[1 3 4 2]);
% 
% % make sure there are only subjects left in the matrix that are actually there and not NaN
% subject_counter = 1;
% for subject = 1:size(all_subjects_prestimulus_baseline_ERSP_dB,4)
%     if ~any(isnan(all_subjects_prestimulus_baseline_ERSP_dB(:,:,:,subject)))
%         all_subjects_prestimulus_baseline_ERSP_dB_noNaN(:,:,:,subject_counter) = all_subjects_prestimulus_baseline_ERSP_dB(:,:,:,subject);
%         all_subjects_ERSP_uncorrected_ERSP_dB_this_cluster_noNaN(:,:,:,subject_counter) = all_subjects_ERSP_uncorrected_ERSP_dB_this_cluster(:,:,:,subject);
%         all_subjects_ERSP_baseline_corrected_dB_noNaN(:,:,:,subject_counter) = all_subjects_ERSP_baseline_corrected_dB_this_cluster(:,:,:,subject);
%         subject_counter = subject_counter+1;
%     end
% end
% 
% % concatenate the data sets per condition, THEN
% % compute significance values  -> the reasoning here is that with null hypothesis uncorrected ERSP
% % data sets should not differ statistically from the baseline data.
% for condition = 1 : length(experiment_conditions_to_test)
%     condition_data_set_for_statistics = {squeeze(all_subjects_ERSP_uncorrected_ERSP_dB_this_cluster_noNaN(condition,:,:,:)) squeeze(all_subjects_prestimulus_baseline_ERSP_dB_noNaN(condition,:,:,:))};
%     [stats, df, p_values(condition,:,:), surrog] = statcond( condition_data_set_for_statistics, 'mode', 'bootstrap', 'naccu', n_permutes);
%     
%     % correct for multiple comparison using false discovery rate
%     [p_values(condition,:,:),p_values_mask(condition,:,:)] = fdr(p_values(condition,:,:),alpha, 'Parametric');
% end
% 
% % statistical difference of the conditions
% difference_data_set_for_statistics = {squeeze(all_subjects_ERSP_baseline_corrected_dB_noNaN(1,:,:,:)) squeeze(all_subjects_ERSP_baseline_corrected_dB_noNaN(2,:,:,:))};
% [stats, df, p_values(3,:,:), surrog] = statcond(difference_data_set_for_statistics, 'mode', 'bootstrap', 'naccu', n_permutes);
% [p_values(3,:,:),p_values_mask(3,:,:)] = fdr(p_values(3,:,:),alpha, 'Parametric');
% 
% 
% disp('statistics calculated')

%     scale_max = max(max(abs(data_to_plot)))/2;
%
%     figure;
%     imagesclogy(times,freqs,data_to_plot,[-scale_max scale_max]);
%     axis xy;
%     cbar;



end