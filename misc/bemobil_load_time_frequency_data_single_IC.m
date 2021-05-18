function time_frequency_data = bemobil_load_time_frequency_data_single_IC(input_path, subject, IC, epochs_info,...
    timewarp_name, trial_normalization, times, timeIndices, freqs, freqIndices, baseline_start_end, experiment_conditions_to_plot,...
    epoch_rejections, epoch_rejections_for_baseline)

filepath_ersp_data = [input_path num2str(subject) '\ERSPs\' timewarp_name '\IC_' num2str(IC) '\'];

load([filepath_ersp_data 'all_epochs_ersp'],'all_epochs_ersp')

all_epochs_ersp = all_epochs_ersp(:, freqIndices(1):freqIndices(2), timeIndices(1):timeIndices(2));

if isempty(epoch_rejections); epoch_rejections = zeros(1,length(epochs_info)); end
if isempty(epoch_rejections_for_baseline); epoch_rejections_for_baseline = zeros(1,length(epochs_info)); end

time_frequency_data.info.baseline_start_end = baseline_start_end;
time_frequency_data.info.epoch_rejections_for_ERSPs = epoch_rejections;
time_frequency_data.info.epoch_rejections_for_baseline = epoch_rejections_for_baseline;
time_frequency_data.info.timewarp_name = timewarp_name;
time_frequency_data.info.experiment_conditions = experiment_conditions_to_plot;
time_frequency_data.freqs = freqs(freqIndices(1):freqIndices(2));
time_frequency_data.times = times(timeIndices(1):timeIndices(2));

% subject_ersp_thisIC_all_epochs_power = NaN(nepochs,length(freqs),length(times));
subject_ersp_thisIC_all_epochs_power = 10.^(all_epochs_ersp/10);
subject_ersp_thisIC_all_epochs_power_unnormalized = subject_ersp_thisIC_all_epochs_power;

if trial_normalization
    full_trial_baselines = nanmean(subject_ersp_thisIC_all_epochs_power,3);
    subject_ersp_thisIC_all_epochs_power = subject_ersp_thisIC_all_epochs_power ./ full_trial_baselines;
end



% assigning a vector with the specified experiment conditions. tests for each epoch if all specified subconditions of a
% condition are met, and if so, assigns that condition number to the epoch. 

if isempty(experiment_conditions_to_plot)
    
    experiment_conditions = ones(1,length(epochs_info));
    
else
    
    experiment_conditions = zeros(1,length(epochs_info));
    
    for n_epoch = 1:length(epochs_info)
        for condition = 1:length(experiment_conditions_to_plot)
            subconditions = false(1,length(experiment_conditions_to_plot(condition).fields));
            
            for subcondition = 1:length(experiment_conditions_to_plot(condition).fields)
                if strcmp(num2str(epochs_info(n_epoch).(experiment_conditions_to_plot(condition).fields{subcondition})),...
                        experiment_conditions_to_plot(condition).conditions{subcondition})
                    
                    subconditions(subcondition) = true;
                    
                end
                
            end
            
            if all(subconditions)
                experiment_conditions(n_epoch) = condition;
            end
            
        end
    end
    
end

% do grand average: if no conditions are specified, all epochs are in the grand average
% if conditions are specified, the grand average should only have the data of the conditions to compare

experiment_conditions_grand_average = logical(experiment_conditions);

% this is the mean ersp of this IC by this subject across all epochs that
% are suitable for the ERSP
subject_ersps_across_epochs_power_grand_average = squeeze(nanmean(subject_ersp_thisIC_all_epochs_power(~logical(epoch_rejections)&experiment_conditions_grand_average,:,:),1));

% this is the mean ersp of this IC by this subject across all epochs that
% are suitable for the baseline of the ERSP
subject_ersps_all_epochs_power_for_baseline_grand_average = subject_ersp_thisIC_all_epochs_power(~logical(epoch_rejections_for_baseline)&experiment_conditions_grand_average,:,:);

% compute baselines of all epochs
subject_baselines_all_epochs_power_for_grand_average = squeeze(nanmean(subject_ersps_all_epochs_power_for_baseline_grand_average(:,:,find(time_frequency_data.times>baseline_start_end(1),1,'first'):find(time_frequency_data.times<=baseline_start_end(2),1,'last')),3));

% compute average baseline
subject_baseline_across_epochs_power_grand_average = squeeze(nanmean(subject_baselines_all_epochs_power_for_grand_average,1))';

time_frequency_data.grand_average.base_power = single(subject_baseline_across_epochs_power_grand_average);
time_frequency_data.grand_average.base_power_ersp = single(repmat(time_frequency_data.grand_average.base_power,1,length(time_frequency_data.times)));
time_frequency_data.grand_average.base_dB_ersp = single(10.*log10(time_frequency_data.grand_average.base_power_ersp)); % this is necessary for bootstrapping
time_frequency_data.grand_average.raw_power = single(subject_ersps_across_epochs_power_grand_average);
time_frequency_data.grand_average.raw_dB = single(10.*log10(subject_ersps_across_epochs_power_grand_average));
time_frequency_data.grand_average.ersp_power = single(subject_ersps_across_epochs_power_grand_average ./ time_frequency_data.grand_average.base_power);
time_frequency_data.grand_average.ersp = single(10.*log10(time_frequency_data.grand_average.ersp_power));
time_frequency_data.grand_average.erspboot = [];
time_frequency_data.grand_average.condition_title = 'grand average';
time_frequency_data.grand_average.n_epochs = sum(~logical(epoch_rejections)&experiment_conditions_grand_average);
time_frequency_data.grand_average.n_epochs_baseline = sum(~logical(epoch_rejections_for_baseline)&experiment_conditions_grand_average);


% Compute unnormalized baseline spectrum for all epochs for later plots. THIS DOES NOT CORRESPOND TO THE ACTUALLY USED
% BASELINES, IF TRIAL NORMALIZATION WAS IN USE!

% this is the mean ersp of this IC by this subject across all epochs that
% are suitable for the baseline of the ERSP
all_epochs_power_for_baseline_grand_average_unnormalized = subject_ersp_thisIC_all_epochs_power_unnormalized(~logical(epoch_rejections_for_baseline)&experiment_conditions_grand_average,:,:);

% compute baselines of all epochs
baselines_all_epochs_power_for_grand_average_unnormalized = squeeze(nanmean(all_epochs_power_for_baseline_grand_average_unnormalized(:,:,find(time_frequency_data.times>baseline_start_end(1),1,'first'):find(time_frequency_data.times<=baseline_start_end(2),1,'last')),3));

time_frequency_data.grand_average.all_epochs_base_power_unnormalized = single(baselines_all_epochs_power_for_grand_average_unnormalized)';

    
for condition = 1:length(experiment_conditions_to_plot)
    
    % this is the mean ersp of this IC by this subject across all epochs that
    % are suitable for the ERSP
    subject_ersps_across_epochs_power_this_condition = squeeze(nanmean(subject_ersp_thisIC_all_epochs_power(~logical(epoch_rejections)&experiment_conditions==condition,:,:,:),1));
    
    % this is the mean ersp of this IC by this subject across all epochs that
    % are suitable for the baseline of the ERSP
    subject_ersps_all_epochs_power_for_baseline_this_condition = subject_ersp_thisIC_all_epochs_power(~logical(epoch_rejections)&experiment_conditions==condition,:,:,:);

    % compute baselines of all epochs
    subject_baselines_all_epochs_power_for_this_condition = squeeze(nanmean(subject_ersps_all_epochs_power_for_baseline_this_condition(:,:,find(time_frequency_data.times>baseline_start_end(1),1,'first'):find(time_frequency_data.times<=baseline_start_end(2),1,'last')),3));

    % compute average baseline
    subject_baseline_across_epochs_power_this_condition = squeeze(nanmean(subject_baselines_all_epochs_power_for_this_condition,1))';

    time_frequency_data.(['condition_' num2str(condition)]).base_power = single(subject_baseline_across_epochs_power_this_condition);
    time_frequency_data.(['condition_' num2str(condition)]).base_power_ersp = single(repmat(time_frequency_data.(['condition_' num2str(condition)]).base_power,1,length(time_frequency_data.times)));
    time_frequency_data.(['condition_' num2str(condition)]).base_dB_ersp = single(10.*log10(time_frequency_data.(['condition_' num2str(condition)]).base_power_ersp)); % this is necessary for bootstrapping
    time_frequency_data.(['condition_' num2str(condition)]).raw_power = single(subject_ersps_across_epochs_power_this_condition);
    time_frequency_data.(['condition_' num2str(condition)]).raw_dB = single(10.*log10(subject_ersps_across_epochs_power_this_condition));
    time_frequency_data.(['condition_' num2str(condition)]).ersp_power = single(subject_ersps_across_epochs_power_this_condition ./ time_frequency_data.(['condition_' num2str(condition)]).base_power);
    time_frequency_data.(['condition_' num2str(condition)]).ersp_with_grand_average_baseline_power = single(subject_ersps_across_epochs_power_this_condition ./ time_frequency_data.grand_average.base_power);
    time_frequency_data.(['condition_' num2str(condition)]).ersp = single(10.*log10(time_frequency_data.(['condition_' num2str(condition)]).ersp_power));
    time_frequency_data.(['condition_' num2str(condition)]).ersp_with_grand_average_baseline = single(10.*log10(time_frequency_data.(['condition_' num2str(condition)]).ersp_with_grand_average_baseline_power));
    
    % Compute unnormalized baseline spectrum for all epochs for later plots. THIS DOES NOT CORRESPOND TO THE ACTUALLY USED
    % BASELINES, IF TRIAL NORMALIZATION WAS IN USE!

    % this is the mean ersp of this IC by this subject across all epochs that
    % are suitable for the baseline of the ERSP
    all_epochs_power_for_baseline_this_cond_unnormalized = subject_ersp_thisIC_all_epochs_power_unnormalized(~logical(epoch_rejections)&experiment_conditions==condition,:,:,:);

    % compute baselines of all epochs
    baselines_all_epochs_power_for_this_cond_unnormalized = squeeze(nanmean(all_epochs_power_for_baseline_this_cond_unnormalized(:,:,find(time_frequency_data.times>baseline_start_end(1),1,'first'):find(time_frequency_data.times<=baseline_start_end(2),1,'last')),3));

    time_frequency_data.(['condition_' num2str(condition)]).all_epochs_base_power_unnormalized = single(baselines_all_epochs_power_for_this_cond_unnormalized)';

    
    % determine title of this condition to save and plot later
    this_condition_title = '';
    this_condition_title = [this_condition_title experiment_conditions_to_plot(condition).fields{1} ': '];
    this_condition_title = [this_condition_title experiment_conditions_to_plot(condition).conditions{1}];
        
    for subcondition = 2:length(experiment_conditions_to_plot(condition).fields)
        this_condition_title = [this_condition_title ', '];
        this_condition_title = [this_condition_title experiment_conditions_to_plot(condition).fields{subcondition} ': '];
        this_condition_title = [this_condition_title experiment_conditions_to_plot(condition).conditions{subcondition}];
        
    end
    
    time_frequency_data.(['condition_' num2str(condition)]).erspboot = [];
    time_frequency_data.(['condition_' num2str(condition)]).condition_title = this_condition_title;
    
    time_frequency_data.(['condition_' num2str(condition)]).n_epochs = sum(~logical(epoch_rejections)&(experiment_conditions==condition));
    time_frequency_data.(['condition_' num2str(condition)]).n_epochs_baseline = sum(~logical(epoch_rejections_for_baseline)&(experiment_conditions==condition));
end
    




% % difference plots
% combinations = nchoosek(experiment_conditions,2);
% combinations_indices = nchoosek(1:length(experiment_conditions),2);
% time_frequency_data.differences.combinations = combinations;
% for combination = 1:size(combinations,1)
% %     time_frequency_data.differences.(['combination_' num2str(combination)]).ersp = single(time_frequency_data.(['condition_' num2str(combinations_indices(combination,1))]).ersp - time_frequency_data.(['condition_' num2str(combinations_indices(combination,2))]).ersp);
% %     time_frequency_data.differences.(['combination_' num2str(combination)]).ersp_power = single(10.^(time_frequency_data.differences.(['combination_' num2str(combination)]).ersp/10)); % power = 10.^(dB/10)
% %     time_frequency_data.differences.(['combination_' num2str(combination)]).base_power = single((time_frequency_data.(['condition_' num2str(combinations_indices(combination,1))]).base_power + time_frequency_data.(['condition_' num2str(combinations_indices(combination,2))]).base_power)/2);
% %     time_frequency_data.differences.(['combination_' num2str(combination)]).base_power_ersp = single(repmat(time_frequency_data.differences.(['combination_' num2str(combination)]).base_power,1,length(times)));
% %     time_frequency_data.differences.(['combination_' num2str(combination)]).base_dB_ersp = single(10.*log10(time_frequency_data.differences.(['combination_' num2str(combination)]).base_power_ersp)); % this is necessary for bootstrapping
% %     time_frequency_data.differences.(['combination_' num2str(combination)]).raw_power = single(time_frequency_data.differences.(['combination_' num2str(combination)]).ersp_power .* time_frequency_data.differences.(['combination_' num2str(combination)]).base_power);
% %     time_frequency_data.differences.(['combination_' num2str(combination)]).raw_dB = single(10.*log10(time_frequency_data.differences.(['combination_' num2str(combination)]).raw_power)); % dB = 10.*log10(power)
% %     time_frequency_data.differences.(['combination_' num2str(combination)]).condition = [combinations{combination,1} ' - ' combinations{combination,2}];
%     
%     time_frequency_data.differences.(['combination_' num2str(combination)]).ersp_power = single(time_frequency_data.(['condition_' num2str(combinations_indices(combination,1))]).ersp_power ./ time_frequency_data.(['condition_' num2str(combinations_indices(combination,2))]).ersp_power);
%     
%     time_frequency_data.differences.(['combination_' num2str(combination)]).ersp = single(10.*log10(time_frequency_data.differences.(['combination_' num2str(combination)]).ersp_power));
%     
%     time_frequency_data.differences.(['combination_' num2str(combination)]).base_power = single(squeeze(mean(time_frequency_data.differences.(['combination_' num2str(combination)]).ersp_power(:,find(times>baseline_start_end(1),1,'first'):find(times<=baseline_start_end(2),1,'last')),2)));
%     time_frequency_data.differences.(['combination_' num2str(combination)]).base_power_ersp = single(repmat(time_frequency_data.differences.(['combination_' num2str(combination)]).base_power,1,length(times)));
%     time_frequency_data.differences.(['combination_' num2str(combination)]).base_dB_ersp = single(10.*log10(time_frequency_data.differences.(['combination_' num2str(combination)]).base_power_ersp)); % this is necessary for bootstrapping
%     time_frequency_data.differences.(['combination_' num2str(combination)]).raw_power = single(time_frequency_data.differences.(['combination_' num2str(combination)]).ersp_power .* time_frequency_data.differences.(['combination_' num2str(combination)]).base_power);
%     time_frequency_data.differences.(['combination_' num2str(combination)]).raw_dB = single(10.*log10(time_frequency_data.differences.(['combination_' num2str(combination)]).raw_power)); % dB = 10.*log10(power)
%     time_frequency_data.differences.(['combination_' num2str(combination)]).condition = [combinations{combination,1} ' - ' combinations{combination,2}];
% end

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