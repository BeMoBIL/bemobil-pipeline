function time_frequency_data = bemobil_load_time_frequency_data_single_IC_single_trials(input_path, subject, IC, epochs_info_filename,...
    timewarp_name, trial_normalization, times, freqs, freqrange, baseline_start_end, baseline_type, latencyMeans,...
    epoch_rejections, epoch_rejections_for_baseline)

% define frequency and time range of ERSP to plot
if isempty(freqrange)
    freqIndices(1) = freqs(1);
    freqIndices(2) = freqs(end);
else
    freqIndices(1) = find(freqs>=freqrange(1),1,'first');
    freqIndices(2) = find(freqs<=freqrange(2),1,'last');
end

% load only with the beginning of the baseline
timeIndices(1) = find(times>=baseline_start_end(1),1,'first');
if isempty(latencyMeans)
    disp('No means in latency present, complete ERSP will be loaded.')
    timeIndices(2) = times(end);
else
    disp('Using last timewarp latency mean as time limit of ERSP')
    timeIndices(2) = find(times<latencyMeans(end),1,'last');
end

% load ersp
filepath_ersp_data = [input_path num2str(subject) '/ERSPs/' timewarp_name '/IC_' num2str(IC) '/'];
load([filepath_ersp_data 'all_epochs_ersp'],'all_epochs_ersp')

% load epoch_info. load stores into a struct, so the first element of the struct has to be taken
epochs_info = load([input_path num2str(subject) '/' epochs_info_filename]);
epoch_info_fields = fieldnames(epochs_info);
epochs_info = epochs_info.(epoch_info_fields{1});

% all_epochs_ersp = all_epochs_ersp(:, freqIndices(1):freqIndices(2), timeIndices(1):timeIndices(2));

if isempty(epoch_rejections); epoch_rejections = zeros(1,length(epochs_info)); end
if isempty(epoch_rejections_for_baseline); epoch_rejections_for_baseline = zeros(1,length(epochs_info)); end

% reject bad trials
ersps_all_epochs_power_for_baseline = all_epochs_ersp(~logical(epoch_rejections_for_baseline),:,:);

% apply baseline
[~, base_corr_ersp_power, ersp_power, ~, ~, base_raw_power, ~] = bemobil_ersp_baseline(...
    ersps_all_epochs_power_for_baseline, times, timeIndices, freqIndices, baseline_start_end, baseline_type, trial_normalization);

% save data & settings
time_frequency_data.uncorrected_ersp_power_single_trials = single(ersp_power);
time_frequency_data.baseline_ersp_power_single_trials = single(base_raw_power);
time_frequency_data.info.baseline_start_end = baseline_start_end;
time_frequency_data.info.epoch_rejections_for_ERSPs = epoch_rejections;
time_frequency_data.info.epoch_rejections_for_baseline = epoch_rejections_for_baseline;
time_frequency_data.info.timewarp_name = timewarp_name;
time_frequency_data.freqs = freqs(freqIndices(1):freqIndices(2));
time_frequency_data.times = times(timeIndices(1):timeIndices(2));

end






