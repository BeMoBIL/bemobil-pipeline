function [average_time_frequency_data, time_frequency_data_of_all_subjects] = bemobil_load_time_frequency_data_of_cluster(STUDY, cluster, input_path, epochs_info_filename,...
    timewarp_name, trial_normalization, times, freqs, baseline_start_end, experiment_conditions_to_plot, freqrange, latencyMeans, n_permutes, alpha, do_movement_based_check, do_auto_epoch_cleaning)



% define frequency and time range of ERSP to plot
if isempty(freqrange)
    freqIndices(1) = freqs(1);
    freqIndices(2) = freqs(end);
else
    freqIndices(1) = find(freqs>freqrange(1),1,'first');
    freqIndices(2) = find(freqs<freqrange(2),1,'last');
end

% load only with the beginning of the baseline
timeIndices(1) = find(times>=baseline_start_end(1),1,'first');
% timeIndices(1) = 1;

if isempty(latencyMeans)
    disp('No means in latency present, complete ERSP will be loaded.')
    
    timeIndices(2) = times(end);
    
else
    disp('Using last timewarp latency mean as time limit of ERSP')
    timeIndices(2) = find(times<latencyMeans(end),1,'last');
end


% the used data sets in the STUDY:
STUDY_sets = cellfun(@str2num, {STUDY.datasetinfo.subject});

unique_setindices = unique(STUDY.cluster(cluster).sets);

unique_subjects = STUDY_sets(unique_setindices);

all_setindices = STUDY.cluster(cluster).sets;

all_sets = STUDY_sets(all_setindices);

all_comps = STUDY.cluster(cluster).comps;

subject_count = 0;

disp(['Loading ERSP data of subjects ' num2str(unique_subjects)])

% load time-freq data per subject
for subject = unique_subjects
    
    subject_count = subject_count + 1;
    
    subject_comps = all_comps(all_sets==subject);
    
    filepath = [input_path num2str(subject) '\'];
    
    % load epoch_info. load stores into a struct, so the first element of the struct has to be taken
    epochs_info = load([filepath '\' epochs_info_filename]);
    epoch_info_fields = fieldnames(epochs_info);
    epochs_info = epochs_info.(epoch_info_fields{1});
    
    
    % THIS IS HARD CODED BULLSHIT! NEXT EXPERIMENT WILL FAIL HERE, SUE ME...
    epoch_rejections_ersp = zeros(1,length(epochs_info));
    epoch_rejections_baseline = zeros(1,length(epochs_info));
    try
        % find out which epochs ought to be rejected based on "wrong" movement during the epoch based on merged checks
        if do_movement_based_check
            epoch_rejections_ersp = for_Marius_extract_removedTrials_MoCap_auto_onsetDetection(subject,epochs_info,[2 3]);
            epoch_rejections_baseline = for_Marius_extract_removedTrials_MoCap_auto_onsetDetection(subject,epochs_info,[1]);
        end
        
    catch
        warning('Epoch rejection based on checks failed! Using all epochs for ERSP and baseline...')
        do_movement_based_check = 0;
    end
    try
        if do_auto_epoch_cleaning
            
            epoch_rejection_auto_cleaning = zeros(1,length(epochs_info));
            
            % edit LG: catch Spot Rotation specific cleaning per cluster
            if contains(filepath, 'Spot_Rotation')
                filename_saveBadEpochIndices=['auto_epoch_cleaning\auto_epoch_cleaning_cluster_'  num2str(cluster)]; % -> THIS MUST MATCH THE BASIC SETTINGS OF THE CLEANING!
                load(['P:\Lukas_Gehrke\studies\Spot_Rotation\data\4_single_subject_analysis\ERSPs\outward\' num2str(subject) '\' filename_saveBadEpochIndices])
                epoch_rejection_auto_cleaning(auto_epoch_cleaning.bad_epochs_final) = 2;
            else
                filename_saveBadEpochIndices= 'epochs_cleaning';
                load([filepath filename_saveBadEpochIndices]);
                epoch_rejection_auto_cleaning(bad_epochs) = 2;
            end
            
            epoch_rejections_ersp = epoch_rejection_auto_cleaning + epoch_rejections_ersp;
            
%             fprintf('Check rejections only: %d\n', sum(epoch_rejections_ersp==1))
%             
%             fprintf('Auto rejections only: %d\n', sum(epoch_rejections_ersp==2))
%             
%             fprintf('Overlapping rejections: %d\n', sum(epoch_rejections_ersp==3))
        end
        
    catch
        warning('Epoch rejection based on automatic cleaning failed! Using all epochs for ERSP and baseline...')
        do_auto_epoch_cleaning = 0;
    end
    
    time_frequency_data_of_all_subjects(subject_count) = bemobil_load_time_frequency_data_single_subject(input_path, subject, subject_comps, epochs_info,...
        timewarp_name, trial_normalization, times, timeIndices, freqs, freqIndices, baseline_start_end, experiment_conditions_to_plot,...
        epoch_rejections_ersp, epoch_rejections_baseline);
    
end

% fill output data

average_time_frequency_data.info = time_frequency_data_of_all_subjects.info;
average_time_frequency_data.info.cluster = cluster;
average_time_frequency_data.info.subjects_used = unique_subjects;
average_time_frequency_data.info.ICs_used = all_comps;
average_time_frequency_data.info.latencyMeans = latencyMeans;
average_time_frequency_data.info.trial_normalization = trial_normalization;
average_time_frequency_data.info.do_auto_epoch_cleaning = do_auto_epoch_cleaning;

average_time_frequency_data.freqs = time_frequency_data_of_all_subjects.freqs;
average_time_frequency_data.times = time_frequency_data_of_all_subjects.times;

% calculate the means

% grand average

% find subjects which have NaN as ERSP, meaning they did not have any epoch of any condition
grand_average_subjec_indices = false(length(time_frequency_data_of_all_subjects),1);
for subject = 1:length(time_frequency_data_of_all_subjects)
    grand_average_subjec_indices(subject) = ~any(any(isnan(time_frequency_data_of_all_subjects(subject).grand_average.ersp)));
end

grand_averages = [time_frequency_data_of_all_subjects(grand_average_subjec_indices).grand_average];

grand_averages_base_power = [grand_averages.base_power];
grand_average_base_power = mean(grand_averages_base_power,2);
average_time_frequency_data.grand_average.base_power = grand_average_base_power;

grand_averages_base_power_ersp = {grand_averages.base_power_ersp};
grand_average_base_power_ersp = mean(cat(3,grand_averages_base_power_ersp{:}),3);
average_time_frequency_data.grand_average.base_power_ersp = grand_average_base_power_ersp;

grand_average_base_dB_ersp = 10.*log10(grand_average_base_power_ersp);
average_time_frequency_data.grand_average.base_dB_ersp = grand_average_base_dB_ersp;

grand_averages_raw_power = {grand_averages.raw_power};
grand_average_raw_power = mean(cat(3,grand_averages_raw_power{:}),3);
average_time_frequency_data.grand_average.raw_power = grand_average_raw_power;

grand_average_raw_dB = 10.*log10(grand_average_raw_power);
average_time_frequency_data.grand_average.raw_dB = grand_average_raw_dB;

grand_average_ersp_power = grand_average_raw_power./grand_average_base_power;
average_time_frequency_data.grand_average.ersp_power = grand_average_ersp_power;

grand_average_ersp = 10.*log10(grand_average_ersp_power);
average_time_frequency_data.grand_average.ersp = grand_average_ersp;

grand_averages_erspboot = [grand_averages.erspboot];
grand_average_erspboot = nanmean(grand_averages_erspboot);
average_time_frequency_data.grand_average.erspboot = grand_average_erspboot;

grand_average_condition = grand_averages(1).condition_title;
average_time_frequency_data.grand_average.condition_title = grand_average_condition;

grand_average_n_epochs = sum([grand_averages.n_epochs]);
average_time_frequency_data.grand_average.n_epochs = grand_average_n_epochs;

grand_average_n_epochs_mean = mean([grand_averages.n_epochs]);
average_time_frequency_data.grand_average.n_epochs_mean = grand_average_n_epochs_mean;

grand_average_n_epochs_std = std([grand_averages.n_epochs]);
average_time_frequency_data.grand_average.n_epochs_std = grand_average_n_epochs_std;

grand_average_n_epochs_baseline_mean = mean([grand_averages.n_epochs_baseline]);
average_time_frequency_data.grand_average.n_epochs_baseline_mean = grand_average_n_epochs_baseline_mean;

grand_average_n_epochs_baseline_std = std([grand_averages.n_epochs_baseline]);
average_time_frequency_data.grand_average.n_epochs_baseline_std = grand_average_n_epochs_baseline_std;

average_time_frequency_data.grand_average.subjects = unique_subjects(grand_average_subjec_indices);

average_time_frequency_data.grand_average.all_epochs_base_power_unnormalized = {grand_averages.all_epochs_base_power_unnormalized};


% this is necessary for bootstrapping
grand_averages_base_dB_all_subjects = 10.*log10(cat(3,grand_averages_base_power_ersp{:}));
grand_averages_raw_dB_all_subjects = 10.*log10(cat(3,grand_averages_raw_power{:}));

% condition(s)
for condition = 1:length(experiment_conditions_to_plot)
    
    % find subjects which have NaN as ERSP, meaning they did not have any epoch of any condition
    this_condition_subjects = false(length(time_frequency_data_of_all_subjects),1);
    for subject = 1:length(time_frequency_data_of_all_subjects)
        this_condition_subjects(subject) = ~any(any(isnan(time_frequency_data_of_all_subjects(subject).(['condition_' num2str(condition)]).ersp)));
    end
    
    this_condition = [time_frequency_data_of_all_subjects(this_condition_subjects).(['condition_' num2str(condition)])];
    
    this_condition_base_power = [this_condition.base_power];
    this_condition_base_power = mean(this_condition_base_power,2);
    average_time_frequency_data.(['condition_' num2str(condition)]).base_power = this_condition_base_power;
    
    this_condition_base_power_ersp_all_sub = {this_condition.base_power_ersp};
    this_condition_base_power_ersp = mean(cat(3,this_condition_base_power_ersp_all_sub{:}),3);
    average_time_frequency_data.(['condition_' num2str(condition)]).base_power_ersp = this_condition_base_power_ersp;
    
    this_condition_base_dB_ersp = 10.*log10(this_condition_base_power_ersp);
    average_time_frequency_data.(['condition_' num2str(condition)]).base_dB_ersp = this_condition_base_dB_ersp;
    
    this_condition_raw_power_all_sub = {this_condition.raw_power};
    this_condition_raw_power = mean(cat(3,this_condition_raw_power_all_sub{:}),3);
    average_time_frequency_data.(['condition_' num2str(condition)]).raw_power = this_condition_raw_power;
    
    this_condition_raw_dB = 10.*log10(this_condition_raw_power);
    average_time_frequency_data.(['condition_' num2str(condition)]).raw_dB = this_condition_raw_dB;
    
    this_condition_ersp_power = this_condition_raw_power./this_condition_base_power;
    average_time_frequency_data.(['condition_' num2str(condition)]).ersp_power = this_condition_ersp_power;
    
    this_condition_ersp = 10.*log10(this_condition_ersp_power);
    average_time_frequency_data.(['condition_' num2str(condition)]).ersp = this_condition_ersp;
    
    this_condition_erspboot = [this_condition.erspboot];
    this_condition_erspboot = nanmean(this_condition_erspboot);
    average_time_frequency_data.(['condition_' num2str(condition)]).erspboot = this_condition_erspboot;
    
    this_condition_condition_title = this_condition(1).condition_title;
    average_time_frequency_data.(['condition_' num2str(condition)]).condition_title = this_condition_condition_title;
    
    this_condition_n_epochs = sum([this_condition.n_epochs]);
    average_time_frequency_data.(['condition_' num2str(condition)]).n_epochs = this_condition_n_epochs;
    
    this_condition_n_epochs_mean = mean([this_condition.n_epochs]);
    average_time_frequency_data.(['condition_' num2str(condition)]).n_epochs_mean = this_condition_n_epochs_mean;
    
    this_condition_n_epochs_std = std([this_condition.n_epochs]);
    average_time_frequency_data.(['condition_' num2str(condition)]).n_epochs_std = this_condition_n_epochs_std;
    
    this_condition_n_epochs_baseline_mean = mean([this_condition.n_epochs_baseline]);
    average_time_frequency_data.(['condition_' num2str(condition)]).n_epochs_baseline_mean = this_condition_n_epochs_baseline_mean;
    
    this_condition_n_epochs_baseline_std = std([this_condition.n_epochs_baseline]);
    average_time_frequency_data.(['condition_' num2str(condition)]).n_epochs_baseline_std = this_condition_n_epochs_baseline_std;
    
    average_time_frequency_data.(['condition_' num2str(condition)]).subjects = unique_subjects(this_condition_subjects);
    
    % Note: The ERSP with grand average baseline (as opposed to the specific condition baseline) can be implemented here
    % as well, but I thought it wasn't necessary...
    
    % this is necessary for bootstrapping
    base_dB_all_subjects.(['condition_' num2str(condition)]) = 10.*log10(cat(3,this_condition_base_power_ersp_all_sub{:}));
    raw_dB_all_subjects.(['condition_' num2str(condition)]) = 10.*log10(cat(3,this_condition_raw_power_all_sub{:}));
    ersp_all_subjects.(['condition_' num2str(condition)]) = 10.*log10(cat(3,this_condition_raw_power_all_sub{:})./cat(3,this_condition_base_power_ersp_all_sub{:}));
    
    average_time_frequency_data.(['condition_' num2str(condition)]).all_epochs_base_power_unnormalized = {this_condition.all_epochs_base_power_unnormalized};

    
end

% difference plot if possible
if length(experiment_conditions_to_plot) == 2
    disp('calculating condition differences...')
    
    difference = struct();
    difference.ersp = average_time_frequency_data.(['condition_' num2str(1)]).ersp - ...
        average_time_frequency_data.(['condition_' num2str(2)]).ersp;
    difference.condition_title = [average_time_frequency_data.(['condition_' num2str(1)]).condition_title ' - ' average_time_frequency_data.(['condition_' num2str(2)]).condition_title];
    
    average_time_frequency_data.difference = difference;
    
end


if ~isempty(n_permutes) && ~n_permutes==0 && ~isempty(alpha)
    disp('Computing significance masks...')
    % grand average
    
    % do permutation test. here the correct freq and time range must be taken since this will affect the values.
    [~, ~, average_time_frequency_data.grand_average.statistics.p_values, ~] =...
        statcond( {grand_averages_base_dB_all_subjects grand_averages_raw_dB_all_subjects}, 'mode', 'bootstrap', 'naccu', n_permutes);
    
    % correct for multiple comparison using false discovery rate
    [average_time_frequency_data.grand_average.statistics.p_value,average_time_frequency_data.grand_average.statistics.p_values_mask] =...
        fdr(average_time_frequency_data.grand_average.statistics.p_values,alpha, 'Parametric');
    
    
    % conditions
    for condition = 1:length(experiment_conditions_to_plot)
        % do permutation test. here the correct freq and time range must be taken since this will affect the values.
        [~, ~, average_time_frequency_data.(['condition_' num2str(condition)]).statistics.p_values, ~] =...
            statcond( {base_dB_all_subjects.(['condition_' num2str(condition)])...
            raw_dB_all_subjects.(['condition_' num2str(condition)])}, 'mode', 'bootstrap', 'naccu', n_permutes);
        
        % correct for multiple comparison using false discovery rate
        [average_time_frequency_data.(['condition_' num2str(condition)]).statistics.p_value,...
            average_time_frequency_data.(['condition_' num2str(condition)]).statistics.p_values_mask] =...
            fdr(average_time_frequency_data.(['condition_' num2str(condition)]).statistics.p_values,alpha,'Parametric');
    end
    
    % difference
    if length(experiment_conditions_to_plot) == 2
        % do permutation test. here the correct freq and time range must be taken since this will affect the values.
        [~, ~, average_time_frequency_data.difference.statistics.p_values, ~] =...
            statcond( {ersp_all_subjects.(['condition_' num2str(1)]) ersp_all_subjects.(['condition_' num2str(2)])},...
            'mode', 'bootstrap', 'naccu', n_permutes);
        
        % correct for multiple comparison using false discovery rate
        [average_time_frequency_data.difference.statistics.p_value,...
            average_time_frequency_data.difference.statistics.p_values_mask] =...
            fdr(average_time_frequency_data.difference.statistics.p_values,alpha, 'Parametric');
    end
    
    % save the significance level
    average_time_frequency_data.alpha = alpha;
    average_time_frequency_data.n_permutes = n_permutes;
    
else
    
    average_time_frequency_data.grand_average.statistics = [];
    
    for condition = 1:length(experiment_conditions_to_plot)
        average_time_frequency_data.(['condition_' num2str(condition)]).statistics = [];
    end
    
    average_time_frequency_data.difference.statistics = [];
end

disp('...done')


