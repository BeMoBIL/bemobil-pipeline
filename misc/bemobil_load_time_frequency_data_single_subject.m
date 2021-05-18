function time_frequency_data = bemobil_load_time_frequency_data_single_subject(input_path, subject, channels, epochs_info,...
    timewarp_name, trial_normalization, times, timeIndices, freqs, freqIndices, baseline_start_end, experiment_conditions_to_plot,...
    epoch_rejections, epoch_rejections_for_baseline,use_channel_data)

if ~exist('use_channel_data','var'); use_channel_data = false; end

channel_count = 0;
fprintf('Subject: %d, channels: ',subject)

for i_channel = channels
    
    channel_count = channel_count + 1;
    
    fprintf('%d, ',i_channel)
    
    if use_channel_data
        time_frequency_data_single_ICs(channel_count) = bemobil_load_time_frequency_data_single_channel(input_path, subject, i_channel,...
            epochs_info, timewarp_name, trial_normalization, times, timeIndices, freqs, freqIndices, baseline_start_end, experiment_conditions_to_plot,...
            epoch_rejections, epoch_rejections_for_baseline);
        
    else
        time_frequency_data_single_ICs(channel_count) = bemobil_load_time_frequency_data_single_IC(input_path, subject, i_channel,...
            epochs_info, timewarp_name, trial_normalization, times, timeIndices, freqs, freqIndices, baseline_start_end, experiment_conditions_to_plot,...
            epoch_rejections, epoch_rejections_for_baseline);
    end
    
end
fprintf('\n')

time_frequency_data.info = time_frequency_data_single_ICs.info;
time_frequency_data.freqs = time_frequency_data_single_ICs.freqs;
time_frequency_data.times = time_frequency_data_single_ICs.times;

% calculate the means

% grand average
grand_averages = [time_frequency_data_single_ICs.grand_average];

grand_averages_base_power = [grand_averages.base_power];
grand_average_base_power = mean(grand_averages_base_power,2);
time_frequency_data.grand_average.base_power = grand_average_base_power;

grand_averages_base_power_ersp = {grand_averages.base_power_ersp};
grand_average_base_power_ersp = mean(cat(3,grand_averages_base_power_ersp{:}),3);
time_frequency_data.grand_average.base_power_ersp = grand_average_base_power_ersp;

grand_average_base_dB_ersp = 10.*log10(grand_average_base_power_ersp); % this is necessary for bootstrapping
time_frequency_data.grand_average.base_dB_ersp = grand_average_base_dB_ersp;

grand_averages_raw_power = {grand_averages.raw_power};
grand_average_raw_power = mean(cat(3,grand_averages_raw_power{:}),3);
time_frequency_data.grand_average.raw_power = grand_average_raw_power;

grand_average_raw_dB = 10.*log10(grand_average_raw_power);
time_frequency_data.grand_average.raw_dB = grand_average_raw_dB;

grand_average_ersp_power = grand_average_raw_power./grand_average_base_power;
time_frequency_data.grand_average.ersp_power = grand_average_ersp_power;

grand_average_ersp = 10.*log10(grand_average_ersp_power);
time_frequency_data.grand_average.ersp = grand_average_ersp;

grand_averages_erspboot = [grand_averages.erspboot];
grand_average_erspboot = nanmean(grand_averages_erspboot);
time_frequency_data.grand_average.erspboot = grand_average_erspboot;

grand_average_condition = grand_averages(1).condition_title;
time_frequency_data.grand_average.condition_title = grand_average_condition;

grand_average_n_epochs = grand_averages(1).n_epochs;
time_frequency_data.grand_average.n_epochs = grand_average_n_epochs;

grand_average_n_epochs_baseline = grand_averages(1).n_epochs_baseline;
time_frequency_data.grand_average.n_epochs_baseline = grand_average_n_epochs_baseline;

time_frequency_data.grand_average.all_epochs_base_power_unnormalized = [grand_averages.all_epochs_base_power_unnormalized];


% conditions
for condition = 1:length(experiment_conditions_to_plot)
    
    this_conditions = [time_frequency_data_single_ICs.(['condition_' num2str(condition)])];
    
    this_conditions_base_power = [this_conditions.base_power];
    this_condition_base_power = mean(this_conditions_base_power,2);
    time_frequency_data.(['condition_' num2str(condition)]).base_power = this_condition_base_power;
    
    this_conditions_base_power_ersp = {this_conditions.base_power_ersp};
    this_condition_base_power_ersp = mean(cat(3,this_conditions_base_power_ersp{:}),3);
    time_frequency_data.(['condition_' num2str(condition)]).base_power_ersp = this_condition_base_power_ersp;
    
    this_condition_base_dB_ersp = 10.*log10(this_condition_base_power_ersp); % this is necessary for bootstrapping
    time_frequency_data.(['condition_' num2str(condition)]).base_dB_ersp = this_condition_base_dB_ersp;
    
    this_conditions_raw_power = {this_conditions.raw_power};
    this_condition_raw_power = mean(cat(3,this_conditions_raw_power{:}),3);
    time_frequency_data.(['condition_' num2str(condition)]).raw_power = this_condition_raw_power;
    
    this_condition_raw_dB = 10.*log10(this_condition_raw_power);
    time_frequency_data.(['condition_' num2str(condition)]).raw_dB = this_condition_raw_dB;
    
    this_condition_ersp_power = this_condition_raw_power./this_condition_base_power;
    time_frequency_data.(['condition_' num2str(condition)]).ersp_power = this_condition_ersp_power;
    
    this_condition_ersp = 10.*log10(this_condition_ersp_power);
    time_frequency_data.(['condition_' num2str(condition)]).ersp = this_condition_ersp;
    
    this_conditions_erspboot = [this_conditions.erspboot];
    this_condition_erspboot = nanmean(this_conditions_erspboot);
    time_frequency_data.(['condition_' num2str(condition)]).erspboot = this_condition_erspboot;
    
    this_condition_condition_title = this_conditions(1).condition_title;
    time_frequency_data.(['condition_' num2str(condition)]).condition_title = this_condition_condition_title;
    
    this_condition_n_epochs = this_conditions(1).n_epochs;
    time_frequency_data.(['condition_' num2str(condition)]).n_epochs = this_condition_n_epochs;
    
    this_condition_n_epochs_baseline = this_conditions(1).n_epochs_baseline;
    time_frequency_data.(['condition_' num2str(condition)]).n_epochs_baseline = this_condition_n_epochs_baseline;
    
    % Note: The ERSP with grand average baseline (as opposed to the specific condition baseline) can be implemented here
    % as well, but I thought it wasn't necessary...
    
    time_frequency_data.(['condition_' num2str(condition)]).all_epochs_base_power_unnormalized = [this_conditions.all_epochs_base_power_unnormalized];
    
end


