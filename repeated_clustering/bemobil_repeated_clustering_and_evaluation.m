% bemobil_repeated_clustering_and_evaluation() - Wrapping function to do the complete repeated clustering and evaluation
% of those clusters on the basis of a region of interest (in MNI coordinates) and a series of weights attributed
% to quality measures of the found ROI clusters. Calls several subfunctions. Saves intermediate data matrices on disk
% (set of all solutions, matrix and plots of multivariate data, top 5 ranked cluster solutions' topoplot and dipoles).
% If the intermediate steps have already been calculated, you can set the flags to do these steps on 0, then the old
% files will be loaded if available. Stores all relevant information in STUDY.etc.bemobil.
%
% Usage:
%   >>  [STUDY, ALLEEG, EEG] = bemobil_repeated_clustering_and_evaluation(STUDY, ALLEEG, EEG, outlier_sigma, n_clust,...
%       n_iterations, cluster_ROI_MNI, quality_measure_weights, do_clustering, do_multivariate_data,...
%       filepath_STUDY, filename_STUDY, filepath_clustering_solutions, filename_clustering_solutions,...
%       filepath_multivariate_data, filename_multivariate_data)
%   
% Inputs:
%   ALLEEG                          - complete EEGLAB data set structure
%   STUDY                           - STUDY data set of EEGLAB, which has to be loaded previously
%   EEG                             - current EEG data set in EEGLAB
%   outlier_sigma                   - standard deviations boundary for outlier detection (e.g. 3)
%   n_clust                         - number of clusters to be created 
%   n_iterations                    - number of iterations to be performed for repeated clustering
%   cluster_ROI_MNI                 - MNI coordinates of the region of interest (THIS NEEDS TO BE A STRUCT consisting of
%                                       .x, .y, and .z fields)
%   quality_measure_weights         - vector of weights for quality measures. 6 entries: 
%                                       subjects, 
%                                       ICs/subjects, 
%                                       normalized spread, 
%                                       mean RV, 
%                                       distance from ROI, 
%                                       mahalanobis distance from median of multivariate distribution (put this very 
%                                           high to get the most "normal" solution)
%   do_clustering                   - [0|1] -> whether or not the clustering should be done (it takes a lot of time). If 0, old
%                                   files will be loaded if possible, if 1, clustering will be redone.
%   do_multivariate_data            - [0|1] -> whether or not the creation of the multivariate dataset should be done. If 0, old
%                                   files will be loaded if possible, if 1, the multivariate dataset will be redone.
%   filepath_STUDY                  - filepath where the new STUDY should be saved
%   filename_STUDY                  - filename of the new STUDY
%   filepath_clustering_solutions   - filepath where the clustering solutions should be saved
%   filename_clustering_solutions   - filename of the repeated clustering solutions' data set
%   filepath_multivariate_data      - filepath where the multivariate data should be saved
%   filename_multivariate_data      - filename of the multivariate data set
%   
%
% Outputs:
%   ALLEEG                          - complete EEGLAB data set structure
%   EEG                             - current EEGLAB EEG structure
%   STUDY                           - STUDY data set of EEGLAB with final best solution as cluster solution
%
%   .study data file of current EEGLAB EEG structure stored on disk 
%   .mat STRUCT data file of the repeated clustering solutions
%   .mat matrix data file of the multivariate data set
%   .fig MATLAB figure of the histograms of the multivariate data
%   .png exports of the figures
%   /evaluations/weights_<quality_measure_weights>/rank-<1-5>_ROI_plot.fig topoplots and dipole plots of the top 5
%       solutions
%   
%
% See also:
%   EEGLAB, bemobil_precluster, bemobil_repeated_clustering, bemobil_create_multivariate_data_from_cluster_solutions, bemobil_clustering_rank_solutions
%
% Authors: Marius Klug, 2021


function [STUDY, ALLEEG, EEG] = bemobil_repeated_clustering_and_evaluation(STUDY, ALLEEG, EEG, outlier_sigma, n_clust,...
    n_iterations, cluster_ROI_MNI, quality_measure_weights, do_clustering, do_multivariate_data, filepath_STUDY,filename_STUDY,...
    filepath_clustering_solutions, filename_clustering_solutions, filepath_multivariate_data, filename_multivariate_data)

% check if files already exist and show warning if it does
% STUDY
mkdir(filepath_STUDY); % make sure that folder exists, nothing happens if so
dir_files = dir(filepath_STUDY);
if ismember(filename_STUDY, {dir_files.name})
    warning([filename_STUDY ' file already exists in: ' filepath_STUDY '. File will be overwritten...']);
end

% clustering_solutions
mkdir(filepath_clustering_solutions); % make sure that folder exists, nothing happens if so
dir_files = dir(filepath_clustering_solutions);

if ismember([filename_clustering_solutions '.mat'], {dir_files.name})
    if do_clustering
        warning([filename_clustering_solutions ' file already exists in: ' filepath_clustering_solutions '. File will be overwritten!!']);
    else
        disp([filename_clustering_solutions ' file already exists in: ' filepath_clustering_solutions]);
    end
else
    disp([filename_clustering_solutions ' does not yet exist in: ' filepath_clustering_solutions '. Repeated clustering will be performed'])
    do_clustering = true;
end

% multivariate_data
mkdir(filepath_multivariate_data); % make sure that folder exists, nothing happens if so
dir_files = dir(filepath_multivariate_data);

if ismember([filename_multivariate_data '.mat'], {dir_files.name})
    if do_multivariate_data
        warning([filename_multivariate_data ' file already exists in: ' filepath_multivariate_data '. File will be overwritten!!']);
    else
        disp([filename_multivariate_data ' file already exists in: ' filepath_multivariate_data]);
    end
else
    disp([filename_multivariate_data ' does not yet exist in: ' filepath_multivariate_data '. Repeated clustering will be performed'])
    do_multivariate_data = true;
end


%% cluster (takes time...)

filename_clustering_solutions_with_path = fullfile(filepath_clustering_solutions, filename_clustering_solutions);

if do_clustering
    
    clustering_solutions = bemobil_repeated_clustering(STUDY,ALLEEG, n_iterations, n_clust, outlier_sigma, STUDY.etc.bemobil.clustering.preclustparams);
    
    disp('Saving clustering solutions...')
    save(filename_clustering_solutions_with_path,'clustering_solutions','-v7.3')
    disp('...done.')
else
    disp('Loading clustering solutions...')
    load(filename_clustering_solutions_with_path,'clustering_solutions')
    disp('...done.')
end

%% create multivariate data of clusters using ROI (takes a little time)

filename_multivariate_data_with_path = fullfile(filepath_clustering_solutions, filename_multivariate_data);

if do_multivariate_data
    
    [cluster_multivariate_data, data_plot] = bemobil_create_multivariate_data_from_cluster_solutions(STUDY,ALLEEG,clustering_solutions,cluster_ROI_MNI);
    
    disp('Saving multivariate data file...')
    save(filename_multivariate_data_with_path,'cluster_multivariate_data','-v7.3')
    set(data_plot,'position',[16 582 1340 751],'color','w')
    saveas(data_plot, [filename_multivariate_data_with_path '_plot.fig'])
    saveas(data_plot, [filename_multivariate_data_with_path '_plot.png'])
    close(data_plot)
    disp('...done.')
    
else
    disp('Loading multivariate data file...')
    load(filename_multivariate_data_with_path,'cluster_multivariate_data')
    disp('...done.')
end

%% find best clustering solutions and store in STUDY. Also create and save top 5 plots

ranked_solutions = bemobil_clustering_rank_solutions(cluster_multivariate_data,quality_measure_weights);

filepath_to_evaluations = fullfile(filepath_clustering_solutions, [filename_multivariate_data '_evaluations'], ['weights_' num2str(quality_measure_weights)]);
mkdir(filepath_to_evaluations)

for rank = 5:-1:1
    
    STUDY.cluster = clustering_solutions.(['solution_' num2str(ranked_solutions(rank))]);
    
    STUDY.etc.bemobil.clustering.parameters = clustering_solutions.parameters;
    STUDY.etc.bemobil.clustering.cluster_ROI_MNI = cluster_multivariate_data.cluster_ROI_MNI;
    STUDY.etc.bemobil.clustering.cluster_ROI_index = cluster_multivariate_data.best_fitting_cluster(ranked_solutions(rank));
    STUDY.etc.bemobil.clustering.cluster_ROI_multivariate_data = cluster_multivariate_data.data(ranked_solutions(rank),:);
    STUDY.etc.bemobil.clustering.cluster_ROI_multivariate_data
    
    topoplot_and_dipplot_fig = figure;
    subplot(1,2,1);
    std_topoplot(STUDY,ALLEEG,'clusters',STUDY.etc.bemobil.clustering.cluster_ROI_index,'figure','off');
    title({get(gca,'title').String
        ['spread = ' num2str(round(STUDY.etc.bemobil.clustering.cluster_ROI_multivariate_data(4),1))]
        ['mean RV = ' num2str(round(STUDY.etc.bemobil.clustering.cluster_ROI_multivariate_data(5)*100,1)) '%']
        ['distance from ROI = ' num2str(round(STUDY.etc.bemobil.clustering.cluster_ROI_multivariate_data(9),1))]})

    subplot(1,2,2);
    
    std_dipplot(STUDY,ALLEEG,'clusters',STUDY.etc.bemobil.clustering.cluster_ROI_index,'figure','off');
    
    set(topoplot_and_dipplot_fig,'position',[16 582 1340 751],'color','w')
    saveas(topoplot_and_dipplot_fig, fullfile(filepath_to_evaluations, ['rank-' num2str(rank) '_ROI_plot.fig']));
    saveas(topoplot_and_dipplot_fig, fullfile(filepath_to_evaluations, ['rank-' num2str(rank) '_ROI_plot.png']));
end

% find dipole locations, centroids, and residual variances of clusters
STUDY = bemobil_dipoles(STUDY,ALLEEG);
STUDY.etc.bemobil.clustering.quality_measure_weights = quality_measure_weights;
fprintf('Best fitting cluster: %d\n', STUDY.etc.bemobil.clustering.cluster_ROI_index)


% save on disk
disp('Saving STUDY...')
[STUDY, EEG] = pop_savestudy( STUDY, EEG, 'filename',filename_STUDY,'filepath',filepath_STUDY);
disp('...done');
eeglab redraw
