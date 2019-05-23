% bemobil_clean_with_iclabel() - Cleans data with the help of the ICLabel classifications. Plots kept and thrown ICs
% according to their respective summed score (or 1 - said score). 
% The classes are: {'Brain'  'Muscle'  'Eye'  'Heart'  'Line Noise'  'Channel Noise'  'Other'}
%
% Usage:
%   >> [ALLEEG, EEG, CURRENTSET] = bemobil_clean_with_iclabel( EEG , ALLEEG, CURRENTSET, classes_to_keep, threshold_to_keep,)
%   >> [ALLEEG, EEG, CURRENTSET] = bemobil_clean_with_iclabel( EEG , ALLEEG, CURRENTSET, classes_to_keep, threshold_to_keep, out_filename, out_filepath)
%
% Inputs:
%   ALLEEG                  - complete EEGLAB data set structure
%   EEG                     - current EEGLAB EEG structure
%   CURRENTSET              - index of current EEGLAB EEG structure within ALLEEG
%   classes_to_keep         - indices of the classes to be kept in the data (e.g. [1 7] for brain and other)
%   threshold_to_keep       - Threshold for the summed score of an IC for it to remain in the dataset (e.g. 0.7)
%   out_filename            - output filename (OPTIONAL ARGUMENT)
%   out_filepath            - output filepath (OPTIONAL ARGUMENT - File will only be saved on disk
%       if both a name and a path are provided)
%
% Outputs:
%   ALLEEG                  - complete EEGLAB data set structure
%   EEG                     - current EEGLAB EEG structure
%   Currentset              - index of current EEGLAB EEG structure within ALLEEG
%
%   .set data file of current EEGLAB EEG structure stored on disk (OPTIONALLY)
%
% See also:
%   EEGLAB, pop_iclabel, bemobil_plot_patterns

function [ALLEEG, EEG, CURRENTSET, ICs_keep, ICs_throw] = bemobil_clean_with_iclabel( EEG , ALLEEG, CURRENTSET, classes_to_keep, threshold_to_keep, out_filename, out_filepath)

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


%% iclabel for artifact rejection

% unnecessary to recompute if present already
if ~isfield(EEG.etc,'ic_classification')

    EEG = pop_iclabel(EEG);

end
% ICs x classes
% classes: {'Brain'  'Muscle'  'Eye'  'Heart'  'Line Noise'  'Channel Noise'  'Other'}

classifications = EEG.etc.ic_classification.ICLabel.classifications;

% brain_and_other = find(any(classifications(:,[1 7])'>0.5));
% artifacts = find(all(classifications(:,[1 7])'<0.5));

summed_scores_to_keep = sum(classifications(:,classes_to_keep),2);
ICs_keep = find(summed_scores_to_keep>threshold_to_keep);

fig1 = bemobil_plot_patterns(EEG.icawinv,EEG.chanlocs,'weights',summed_scores_to_keep,'minweight',threshold_to_keep);
% title(['Classes to keep: ' num2str(classes_to_keep) ', Threshold: ' num2str(threshold_to_keep)]); % doesn't work since it titles only one the last plot, not the figure


ICs_throw = find(summed_scores_to_keep<threshold_to_keep);
summed_scores_to_throw = 1 - summed_scores_to_keep;

fig2 = bemobil_plot_patterns(EEG.icawinv,EEG.chanlocs,'weights',summed_scores_to_throw,'minweight',1 - threshold_to_keep);
% title(['Classes to throw out: ' num2str(artifacts_from_summed) ', Threshold: ' num2str(1 - threshold_to_keep)]);


EEG_clean = pop_subcomp( EEG, ICs_throw);

EEG_clean.etc.ic_cleaning.classes_to_keep = classes_to_keep;
EEG_clean.etc.ic_cleaning.threshold_to_keep = threshold_to_keep;
EEG_clean.etc.ic_cleaning.summed_scores_to_keep = summed_scores_to_keep;
EEG_clean.etc.ic_cleaning.ICs_keep = ICs_keep;
EEG_clean.etc.ic_cleaning.ICs_throw = ICs_throw;


%% new data set in EEGLAB
[ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG_clean, CURRENTSET, 'gui', 'off');
EEG = eeg_checkset( EEG );

% save on disk
if save_file_on_disk
    EEG = pop_saveset( EEG, 'filename',out_filename,'filepath', out_filepath);
    disp('...done');
	savefig(fig1,fullfile(out_filepath,'ICs_kept'))
	savefig(fig2,fullfile(out_filepath,'ICs_thrown_out'))
end