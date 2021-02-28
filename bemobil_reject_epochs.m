function EEG_epoched = bemobil_reject_epochs(EEG_epoched,fixed_threshold,weights,...
	use_max_epochs,max_epochs,use_kneepoint,kneepoint_offset,do_plot)
%
% color scheme was checked with https://www.color-blindness.com/coblis-color-blindness-simulator/
%
% SEE ALSO: BEMOBIL_REJECT_CONTINUOUS

%% features
% simple, fast, but surprisingly informative... best to remove
% eye-components before, or to filter with a cutoff (not passband!) > 12 Hz

means = squeeze(nanmean(EEG_epoched.data,2));
SDs = squeeze(nanstd(EEG_epoched.data,[],2));

mean_of_means = mean(abs(means));
SDs_of_means = std(means);
mean_of_SDs = mean(SDs);
SDs_of_SDs = std(SDs);

all_normalized_features = [...
	mean_of_means/median(mean_of_means)
	SDs_of_means/median(SDs_of_means)
	mean_of_SDs/median(mean_of_SDs)
	SDs_of_SDs/median(SDs_of_SDs)];

weighted_sum_features = sum(all_normalized_features .* weights(:));

%% determine rejection based on weighted features

sorted = sort(weighted_sum_features);

if use_kneepoint
	% abuse an image processing function i found on the interwebz to find the
	% knee point here. works better than the knee_pt function imo.
	% https://de.mathworks.com/matlabcentral/answers/483969-find-knee-elbow-of-a-curve
	% author is Mark Hayworth!
% 	kneepoint = triangle_threshold(sorted, 'L', 0);
% 	rejection_index = kneepoint - round(length(sorted)/100*kneepoint_offset);
% 	reject_epochs = weighted_sum_features > sorted(rejection_index);
    
	[n_remove, threshold] = bemobil_iterative_threshold_detection(weighted_sum_features);
    kneepoint = length(weighted_sum_features)-n_remove;
	rejection_index = kneepoint - round(length(sorted)/100*kneepoint_offset);
    
elseif use_max_epochs
	rejection_index = max_epochs;
else
	rejection_index = round(length(weighted_sum_features)*(1-fixed_threshold));
	
	% median based rejection
% 	fixed_threshold = median(sorted)*2-sorted(1);
% 	reject_epochs = weighted_sum_features > fixed_threshold;
	
end

	reject_epochs = weighted_sum_features > sorted(rejection_index);
	
%% plot

if do_plot
	
	figure; set(gcf,'color','w','position',[0 0 1920 1080])
	
	% dataset with 11 channels (to not clutter the plot and make it
	% unnecessarily large)
	subplot(2,5,[1:5]); hold on
	selected_channels = round(linspace(1,EEG_epoched.nbchan,11));
	all_data_connected = reshape(EEG_epoched.data,EEG_epoched.nbchan,EEG_epoched.pnts*EEG_epoched.trials);
	all_samples = 1:EEG_epoched.pnts*EEG_epoched.trials;
	
	% plot all epochs in blue
	plot(all_samples/EEG_epoched.srate,all_data_connected(selected_channels,:)','color',[0 0.5 0.8])

	% overwrite rejected epochs in red - this is actually more efficient
	% than plotting all epochs either in blue or red, just MATLAB things i
	% guess...
	for n = 1:size(EEG_epoched.data,3)
		
		lower_bound = 1+(n-1)*size(EEG_epoched.data,2);
		upper_bound = min(EEG_epoched.pnts*EEG_epoched.trials, lower_bound+size(EEG_epoched.data,2)-1);

		if reject_epochs(n)
			plot(all_samples(lower_bound:upper_bound)/EEG_epoched.srate,...
				EEG_epoched.data(selected_channels,lower_bound:upper_bound)','color',[0.7 0 0])
		end
		
	end
	
	title(['Concatenated Epochs (not continuous data!), rejected = '...
		num2str(round(sum(reject_epochs)/length(reject_epochs)*100,2)) '%'])
	xlim([0 length(all_samples)/EEG_epoched.srate])
	
	
	% Weighted Sum
	subplot(256); hold on
	
	sorted = sort(weighted_sum_features);
	plot(sorted,'linewidth',2);
	xlim([0 length(reject_epochs)])
	ylim([0 max(sorted)])
	plot([rejection_index rejection_index],[0 max(sorted)],'color',[0.7 0 0],'linewidth',2);
	if use_kneepoint
		plot(kneepoint,sorted(kneepoint),'k*','markersize',15);
		title({['Sums of weighted normalized scores']
			['kneepoint offset = ' num2str(kneepoint_offset) '%']})
	else
		title({['Sums of weighted normalized scores']
			['fixed threshold = ' num2str(fixed_threshold * 100) '%']})
	end
	axis square
	grid on
	
	% Mean of Means
	subplot(257); hold on
	
	sorted = sort(mean_of_means);
	plot(sorted,'linewidth',2);
	xlim([0 length(reject_epochs)])
	ylim([0 max(sorted)])
	plot([rejection_index rejection_index],[0 max(sorted)],'color',[0.7 0 0],'linewidth',2);
	title(['Mean of Means, weight = ' num2str(weights(1))])
	axis square
	grid on
	
	% SDs of Means
	subplot(258); hold on
	
	sorted = sort(SDs_of_means);
	plot(sorted,'linewidth',2);
	xlim([0 length(reject_epochs)])
	ylim([0 max(sorted)])
	plot([rejection_index rejection_index],[0 max(sorted)],'color',[0.7 0 0],'linewidth',2);
	title(['SDs of Means, weight = ' num2str(weights(2))])
	axis square
	grid on
	
	% Mean of SDs
	subplot(259); hold on
	
	sorted = sort(mean_of_SDs);
	plot(sorted,'linewidth',2);
	xlim([0 length(reject_epochs)])
	ylim([0 max(sorted)])
	plot([rejection_index rejection_index],[0 max(sorted)],'color',[0.7 0 0],'linewidth',2);
	title(['Mean of SDs, weight = ' num2str(weights(3))])
	axis square
	grid on
	
	% SDs of SDs
	subplot(2,5,10); hold on
	
	sorted = sort(SDs_of_SDs);
	plot(sorted,'linewidth',2);
	xlim([0 length(reject_epochs)])
	ylim([0 max(sorted)])
	plot([rejection_index rejection_index],[0 max(sorted)],'color',[0.7 0 0],'linewidth',2);
	title(['SDs of SDs, weight = ' num2str(weights(4))])
	axis square
	grid on
	
end

%% store info 

EEG_epoched.etc.bemobil_reject_epochs.rejected_epochs = reject_epochs;
EEG_epoched.etc.bemobil_reject_epochs.proportion_rejected = round(sum(reject_epochs)/length(reject_epochs)*100,2);
EEG_epoched.etc.bemobil_reject_epochs.fixed_threshold = fixed_threshold;
EEG_epoched.etc.bemobil_reject_epochs.weights = weights;
EEG_epoched.etc.bemobil_reject_epochs.use_max_epochs = use_max_epochs;
EEG_epoched.etc.bemobil_reject_epochs.max_epochs = max_epochs;
EEG_epoched.etc.bemobil_reject_epochs.use_kneepoint = use_kneepoint;
EEG_epoched.etc.bemobil_reject_epochs.kneepoint_offset = kneepoint_offset;
EEG_epoched.etc.bemobil_reject_epochs.weighted_sum_features = weighted_sum_features;
EEG_epoched.etc.bemobil_reject_epochs.mean_of_means = mean_of_means;
EEG_epoched.etc.bemobil_reject_epochs.SDs_of_means = SDs_of_means;
EEG_epoched.etc.bemobil_reject_epochs.mean_of_SDs = mean_of_SDs;
EEG_epoched.etc.bemobil_reject_epochs.SDs_of_SDs = SDs_of_SDs;
