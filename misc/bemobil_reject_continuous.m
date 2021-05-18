% BEMOBIL_REJECT_CONTINUOUS - Reject continuous data automatically.
% 
% The data is being split in epochs which are then automatically rejected.
% Using mean of channel means (to detect general large shifts), mean of channel SDs (to
% detect large fluctuations in many channels, e.g. large muscle activity),
% SD of channel means (to detect phases where single channels went off), and SD of channel SDs
% (to detect large fluctuations in single channels, even when not going off
% in one direction) as features. 
% 
% All info is stored in EEG.etc.bemobil_reject_continuous.
%
% Usage:
%	[ ALLEEG EEG_cleaned CURRENTSET ] = bemobil_reject_continuous(ALLEEG, EEG, CURRENTSET,...
% 	epochs_length, epochs_overlap, epoch_buffer, fixed_threshold, weights,...
% 	use_kneepoint, kneepoint_offset, highpass_cutoff, do_plot)
% 
% Inputs:
%   ALLEEG             - complete EEGLAB data set structure
%   EEG                - current EEGLAB EEG structure
%   CURRENTSET         - index of current EEGLAB EEG structure within ALLEEG
%   epochs_length      - length of the generated epochs (in s, e.g. 0.5)
%   epochs_overlap     - overlap of the generated epochs (in s, e.g. 0.125)
%   epoch_buffer       - buffer which is applied around rejected epochs (in s, e.g. 0.0625)
%   fixed_threshold    - threshold in proportion of epochs to reject (e.g. 0.07)
%   weights            - weights for the different features. recommended to keep at [1 1 1 1]
%   use_kneepoint      - if 1, automatic kneepoint detection is used
%								instead of fixed threshold for rejection.
%								can lead to issues if data is too clean or
%								too noisy because the knee point is shifted
%								and does not make sense any more. not
%								recommended, just here for testing and
%								demonstration purposes.
%   kneepoint_offset   - offset in % points added to the rejected knee
%								point threshold to make it more stringent
%								or lax.
%   highpass_cutoff    - filter frequency to use before detecting
%								artifacts cutoff frequency in Hz, e.g. 10). 
%								it is recommended to set this to a value
%								that removes eye artifacts, since they can
%								be removed by ICA and should not be removed
%								in the time domain.
%   do_plot            - plot diagnostics figures for epoch cleaning
%	do_plot_continuous - visualize continuous artifact rejection
% Outputs:
%   ALLEEG                  - complete EEGLAB data set structure
%   EEG_cleaned             - cleaned EEG dataset
%   Currentset              - index of current EEGLAB EEG structure within ALLEEG
%
% See also:
%   eeglab, bemobil_reject_epochs
%
% Authors: Marius Klug, 2020

function [ ALLEEG, EEG_cleaned, CURRENTSET, plot_handles ] = bemobil_reject_continuous(ALLEEG, EEG, CURRENTSET,...
	epochs_length, epochs_overlap, epoch_buffer, fixed_threshold, weights,...
	use_kneepoint, kneepoint_offset, highpass_cutoff, do_plot, do_plot_continuous)

%% input checks 
% (dirty for now because i want to change the entire system)

if ~ismatrix(EEG.data)
	error('Cannot call this funtion on an epoched dataset, call bemobil_reject_epochs directly instead!')
end

if ~exist('do_plot_continuous')
	do_plot_continuous = 1;
end

if isempty(highpass_cutoff) 
	
	% default filter cutoff to remove eye artifacts in order to disregard
	% them during cleaning
	highpass_cutoff = 10;
	
end

%% filter before epoching, if no bemobil filter info was found
EEG_for_cleaning = EEG;

% remove events to save computation and RAM
EEG_for_cleaning.event = [];
EEG_for_cleaning.urevent = [];
	
if highpass_cutoff == 0
	
	disp('Filtering is disabled ("filter_freq" set to 0).')
	
elseif ~isfield(EEG.etc, 'filter') || ~isfield(EEG.etc.filter.highpass, 'cutoff')
	
	disp(['No "cutoff" filter information was found in "EEG.etc.filter.highpass".'])
	disp(['Filtering before automatic cleaning. This will not be applied to the final dataset!'])
	disp('This can be disabled by setting the "filter_freq" to 0.')
	% with an order of 414 the transition bandwidth is 2 Hz, so the
	% cutoff+1 is the passband
	[ ALLEEG EEG_for_cleaning CURRENTSET ] = bemobil_filter(ALLEEG, EEG_for_cleaning, CURRENTSET, highpass_cutoff+1, [],[],[],412);
	
elseif EEG.etc.filter.highpass.cutoff < highpass_cutoff
	
	disp(['The filter cutoff that was applied previously (' num2str(EEG.etc.filter.highpass.cutoff)...
		' Hz, found in "EEG.etc.filter.highpass") was not as high as the defined "filter_freq" (' num2str(highpass_cutoff) ' Hz).'])
	disp(['Filtering again before automatic cleaning. This will not be applied to the final dataset!'])
	disp('This can be disabled by setting the "filter_freq" to 0.')
	[ ALLEEG EEG_for_cleaning CURRENTSET ] = bemobil_filter(ALLEEG, EEG_for_cleaning, CURRENTSET, highpass_cutoff+1, [],[],[],412);
	
else
	
	disp('Data already sufficiently filtered.')
	
end


%% create epochs

[EEG_for_cleaning] = eeg_regepochs(EEG_for_cleaning, epochs_length-epochs_overlap,...
	[0 epochs_length], NaN);


%% call rejection algorithm on epoched dataset

EEG_for_cleaning = bemobil_reject_epochs(EEG_for_cleaning,fixed_threshold,weights,...
	0,0,use_kneepoint,kneepoint_offset,do_plot);

if do_plot
	plot_handles(1) = gcf;
end

%% find original EEG indices

clean_sample_mask = ones(1,EEG.pnts);

for n = 1:size(EEG_for_cleaning.data,3)

	urevent = EEG_for_cleaning.epoch(n).eventurevent{1};
	
	% latency is given in samples, not ms!
	lower_bound = round(max(1,EEG_for_cleaning.urevent(urevent).latency-EEG.srate*epoch_buffer));
	
	% if it's the last epoch, continue to end of dataset, otherwise epoch
	% length
	upper_bound = round(min(EEG.pnts,lower_bound+size(EEG_for_cleaning.data,2)-1+2*EEG_for_cleaning.srate*epoch_buffer));

	clear clean_vector_this_epoch
	clean_vector_this_epoch(1:upper_bound-lower_bound+1) = ...
		deal(~EEG_for_cleaning.etc.bemobil_reject_epochs.rejected_epochs(n));
	
	clean_sample_mask(lower_bound:upper_bound) =...
		clean_vector_this_epoch & clean_sample_mask(lower_bound:upper_bound);
end

clean_sample_mask = logical(clean_sample_mask);

%% find latency of regions and reject
retain_data_intervals = reshape(find(diff([false clean_sample_mask false])),2,[])';
retain_data_intervals(:,2) = retain_data_intervals(:,2)-1;

% reject regions from original dataset
EEG_cleaned = pop_select(EEG, 'point', retain_data_intervals);
EEG_cleaned.etc.clean_sample_mask = clean_sample_mask;

%% visualize
if do_plot_continuous
	vis_artifacts(EEG_cleaned,EEG);
	plot_handles(2) = gcf;
end

%% store info
EEG_cleaned.etc.bemobil_reject_continuous.proportion_rejected_samples = (1 - sum(clean_sample_mask) / length(clean_sample_mask))*100;
EEG_cleaned.etc.bemobil_reject_continuous.clean_samples = clean_sample_mask;
EEG_cleaned.etc.bemobil_reject_continuous.retain_data_intervals = retain_data_intervals;
EEG_cleaned.etc.bemobil_reject_continuous.applied_filter_before_cleaning = EEG_for_cleaning.etc.filter;
EEG_cleaned.etc.bemobil_reject_continuous.epochs_length = epochs_length;
EEG_cleaned.etc.bemobil_reject_continuous.epochs_overlap = epochs_overlap;
EEG_cleaned.etc.bemobil_reject_continuous.epoch_buffer = epoch_buffer;
EEG_cleaned.etc.bemobil_reject_continuous.fixed_threshold = fixed_threshold;
EEG_cleaned.etc.bemobil_reject_continuous.weights = weights;
EEG_cleaned.etc.bemobil_reject_continuous.use_kneepoint = use_kneepoint;
EEG_cleaned.etc.bemobil_reject_continuous.kneepoint_offset = kneepoint_offset;
EEG_cleaned.etc.bemobil_reject_continuous.bemobil_reject_epochs = EEG_for_cleaning.etc.bemobil_reject_epochs;


