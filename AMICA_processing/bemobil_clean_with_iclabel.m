% bemobil_clean_with_iclabel() - Cleans data with the help of the ICLabel classifications. Plots ICs that are kept
% scaled to their respective summed score.
% The classes are: {'Brain'  'Muscle'  'Eye'  'Heart'  'Line Noise'  'Channel Noise'  'Other'}
%
% Usage:
%   >> [ALLEEG, EEG, CURRENTSET, ICs_keep, ICs_throw] = bemobil_clean_with_iclabel( EEG , ALLEEG, CURRENTSET, classifier_version,...
%     classes_to_keep, threshold_to_keep, out_filename, out_filepath)
% Inputs:
%   ALLEEG                  - complete EEGLAB data set structure
%   EEG                     - current EEGLAB EEG structure
%   CURRENTSET              - index of current EEGLAB EEG structure within ALLEEG
%   classifier_version      - ICLabel classifier version to use (DEFAULT = 'lite')
%   classes_to_keep         - indices of the classes to be kept in the data (e.g. [1 7] for brain and other, DEFAULT = 1)
%   threshold_to_keep       - Threshold for the summed score of an IC for it to remain in the dataset (e.g. 0.7). If -1,
%                               popularity vote will be used (DEFAULT = -1)
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
% 
% Author: Marius Klug, 2021

function [ALLEEG, EEG, CURRENTSET, ICs_keep, ICs_throw] = bemobil_clean_with_iclabel( EEG , ALLEEG, CURRENTSET, classifier_version,...
    classes_to_keep, threshold_to_keep, out_filename, out_filepath)

if ~exist('classifier_version','var') || isempty(classifier_version)
    classifier_version = 'lite';
end

if ~exist('classes_to_keep','var') || isempty(classes_to_keep)
    classes_to_keep = 1;
end

if ~exist('threshold_to_keep','var') || isempty(threshold_to_keep)
    threshold_to_keep = -1;
end

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

% add to path for plotting
addpath(fullfile(fileparts(which('iclabel')),'viewprops'))

%% iclabel for artifact rejection

% unnecessary to recompute if present already
if ~isfield(EEG.etc,'ic_classification') || ~strcmp(EEG.etc.ic_classification.ICLabel.version,classifier_version)
    
    EEG = iclabel(EEG,classifier_version);
    
end
% ICs x classes
% classes: {'Brain'  'Muscle'  'Eye'  'Heart'  'Line Noise'  'Channel Noise'  'Other'}

classifications = EEG.etc.ic_classification.ICLabel.classifications;
[classes ~] = find(transpose(EEG.etc.ic_classification.ICLabel.classifications == max(EEG.etc.ic_classification.ICLabel.classifications,[],2)));

% brain_and_other = find(any(classifications(:,[1 7])'>0.5));
% artifacts = find(all(classifications(:,[1 7])'<0.5));


if threshold_to_keep ~= -1
    summed_scores_to_keep = sum(classifications(:,classes_to_keep),2);
    ICs_keep = find(summed_scores_to_keep>threshold_to_keep);
    ICs_throw = find(summed_scores_to_keep<threshold_to_keep);
    
    summed_scores_to_keep(ICs_throw) = 0;
%     pop_viewprops(EEG, 0,ICs_keep');
%     fig1 = gcf;
%     pop_viewprops(EEG, 0,ICs_throw');
%     fig2 = gcf;

    titles = {};
    for i_title = 1:size(EEG.icawinv,2)
        titles{i_title} = [num2str(i_title) ': ' num2str(round(summed_scores_to_keep(i_title),3))];
    end
    fig1 = bemobil_plot_patterns(EEG.icawinv,EEG.chanlocs,'weights',summed_scores_to_keep,'minweight',0.01,'titles',titles);
    
%     pop_viewprops(EEG, 0,ICs_keep');
%     fig1 = gcf;
%     fig1 = bemobil_plot_patterns(EEG.icawinv,EEG.chanlocs,'weights',summed_scores_to_keep,'minweight',threshold_to_keep);
%     % title(['Classes to keep: ' num2str(classes_to_keep) ', Threshold: ' num2str(threshold_to_keep)]); % doesn't work since it titles only one the last plot, not the figure
%     
%     fig2 = bemobil_plot_patterns(EEG.icawinv,EEG.chanlocs,'weights',summed_scores_to_throw,'minweight',1 - threshold_to_keep);
%     % title(['Classes to throw out: ' num2str(artifacts_from_summed) ', Threshold: ' num2str(1 - threshold_to_keep)]);
else
    ICs_keep = zeros(length(classes),1);
    summed_scores_to_keep = zeros(size(EEG.icawinv,2),1);
    for class_to_keep = classes_to_keep
        ICs_keep(classes == class_to_keep) = 1;
        summed_scores_to_keep = summed_scores_to_keep + EEG.etc.ic_classification.ICLabel.classifications(:,1);
    end
    
    ICs_throw = find(~ICs_keep);
    ICs_keep = find(ICs_keep);
    
    summed_scores_to_keep(ICs_throw) = 0;
%     pop_viewprops(EEG, 0,ICs_keep');
%     fig1 = gcf;
%     pop_viewprops(EEG, 0,ICs_throw');
%     fig2 = gcf;

    titles = {};
    for i_title = 1:size(EEG.icawinv,2)
        titles{i_title} = [num2str(i_title) ': ' num2str(round(summed_scores_to_keep(i_title),3))];
    end
    fig1 = bemobil_plot_patterns(EEG.icawinv,EEG.chanlocs,'weights',summed_scores_to_keep,'minweight',0.01,'titles',titles);
    % title(['Classes to keep: ' num2str(classes_to_keep) ', Threshold: ' num2str(threshold_to_keep)]); % doesn't work since it titles only one the last plot, not the figure
    
%     EEG_clean = pop_subcomp( EEG, ICs_keep);
%     fig2 = bemobil_plot_patterns(EEG_clean.icawinv,EEG_clean.chanlocs,'weights',summed_scores_to_throw(ICs_throw));
    % title(['Classes to throw out: ' num2str(artifacts_from_summed) ', Threshold: ' num2str(1 - threshold_to_keep)]);
end

EEG_clean = pop_subcomp( EEG, ICs_throw);
set(fig1,'color','w','position',get(0,'screensize'))

% plot dipoles
pop_dipplot( EEG_clean, 1:size(EEG_clean.icaweights,1),...
    'mri',fullfile(fileparts(which('dipfitdefs')), 'standard_BEM','standard_mri.mat'),'normlen','on');
view(30,30)

%%
EEG_clean.etc.ic_cleaning.classes_to_keep = classes_to_keep;
EEG_clean.etc.ic_cleaning.threshold_to_keep = threshold_to_keep;
EEG_clean.etc.ic_cleaning.ICs_keep = ICs_keep;
EEG_clean.etc.ic_cleaning.ICs_throw = ICs_throw;

%% new data set in EEGLAB
[ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG_clean, CURRENTSET, 'gui', 'off');
EEG = eeg_checkset( EEG );

% save on disk
if save_file_on_disk
    
    % save dipole fig
    print(gcf,fullfile(out_filepath,[erase(out_filename,'.set') '_brain_dipoles.png']),'-dpng')
    savefig(fullfile(out_filepath,[erase(out_filename,'.set') '_brain_dipoles.fig']))

    close
    print(fig1,fullfile(out_filepath,[erase(out_filename,'.set') '_ICs_kept.png']),'-dpng')
    close
    EEG = pop_saveset( EEG, 'filename',erase(out_filename,'.set'),'filepath', out_filepath);
    disp('...done');
end