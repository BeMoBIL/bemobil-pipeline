function [ALLEEG, EEG, CURRENTSET] = bemobil_ssd_frontal_parietal(ALLEEG, EEG, CURRENTSET, bemobil_config, out_filename, out_filepath)

% only save a file on disk if both a name and a path are provided
save_file_on_disk = (exist('out_filename', 'var') && exist('out_filepath', 'var'));

% check if file already exist and show warning if it does
if save_file_on_disk
	mkdir(out_filepath); % make sure that folder exists, nothing happens if so
	dir_files = dir(out_filepath);
	if ismember(out_filename, {dir_files.name})
		warning([out_filename ' file already exists in: ' out_filepath '. File will be overwritten...']);
	end
end

%% find channel indices

frontal_chanindices = [];

for i_frontal_channame = 1:length(bemobil_config.frontal_channames)
	frontal_chanindices(end+1) = find(strcmp({EEG.chanlocs.labels},bemobil_config.frontal_channames{i_frontal_channame}));
end

parietal_chanindices = [];

for i_parietal_channame = 1:length(bemobil_config.parietal_channames)
	parietal_chanindices(end+1) = find(strcmp({EEG.chanlocs.labels},bemobil_config.parietal_channames{i_parietal_channame}));
end

%% region 1 (frontal theta)

[ALLEEG EEG_theta CURRENTSET] = bemobil_signal_decomposition_extended(ALLEEG, EEG, CURRENTSET,...
	'SSD', 'ssd_freq', bemobil_config.ssd_freq_theta);

frontal_weights = EEG_theta.icawinv(frontal_chanindices,:);

scaled_frontal_weights = abs(frontal_weights ./ max(abs(EEG_theta.icawinv),[],1));

parietal_weights = EEG_theta.icawinv(parietal_chanindices,:);

scaled_parietal_weights = abs(parietal_weights ./ max(abs(EEG_theta.icawinv),[],1));

frontal_theta_pattern = find(sum(scaled_frontal_weights,1) == max(sum(scaled_frontal_weights - scaled_parietal_weights,1)));

% pop_topoplot(EEG_theta,0, [1:15] ,'Merged datasets resampled',[3 5] ,0,'electrodes','on');
% pop_prop( EEG_theta, 0, frontal_theta_pattern, NaN, {'freqrange' [2 50] });

%% region 2 (parietal alpha)

[ALLEEG EEG_alpha CURRENTSET] = bemobil_signal_decomposition_extended(ALLEEG, EEG, CURRENTSET,...
	'SSD', 'ssd_freq', bemobil_config.ssd_freq_alpha);

parietal_weights = EEG_alpha.icawinv(parietal_chanindices,:);

scaled_parietal_weights = abs(parietal_weights ./ max(abs(EEG_alpha.icawinv),[],1));

frontal_weights = EEG_alpha.icawinv(frontal_chanindices,:);

scaled_frontal_weights = abs(frontal_weights ./ max(abs(EEG_alpha.icawinv),[],1));

parietal_alpha_pattern = find(sum(scaled_parietal_weights-scaled_frontal_weights,1) == max(sum(scaled_parietal_weights-scaled_frontal_weights,1)));

% pop_topoplot(EEG_alpha,0, [1:15] ,'Merged datasets resampled',[3 5] ,0,'electrodes','on');
% pop_prop( EEG_alpha, 0, parietal_alpha_pattern, NaN, {'freqrange' [2 50] });
figure; plot(sum(scaled_parietal_weights,1))

%% store only these two components
EEG.icawinv(:,1) = EEG_theta.icawinv(:,frontal_theta_pattern);
EEG.icawinv(:,2) = EEG_alpha.icawinv(:,parietal_alpha_pattern);
EEG.icaweights = eye(2);
EEG.icasphere(1,:) = EEG_theta.icasphere(frontal_theta_pattern,:);
EEG.icasphere(2,:) = EEG_alpha.icasphere(parietal_alpha_pattern,:);
EEG.icachansind = 1:EEG.nbchan;

EEG.etc.spatial_filter.algorithm = 'SSD_frontal_parietal';
EEG.etc.spatial_filter.ssd.frontal_channames = bemobil_config.frontal_channames;
EEG.etc.spatial_filter.ssd.parietal_channames = bemobil_config.parietal_channames;
EEG.etc.spatial_filter.ssd.ssd_freq_frontal = bemobil_config.ssd_freq_theta;
EEG.etc.spatial_filter.ssd.ssd_freq_parietal = bemobil_config.ssd_freq_alpha;

EEG = eeg_checkset(EEG);

% new data set in EEGLAB
[ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, CURRENTSET, 'gui', 'off');
EEG = eeg_checkset( EEG );

% save on disk
if save_file_on_disk
	EEG = pop_saveset( EEG, 'filename',out_filename,'filepath', out_filepath);
	disp('...done');
end

[ALLEEG EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);
