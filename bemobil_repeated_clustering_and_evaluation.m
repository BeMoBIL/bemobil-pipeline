% bemobil_repeated_clustering_and_evaluation() - Wrapping function to do the complete repeated clustering and evaluation
% of those clusters on the basis of a region of interest (in talairach coordinates) and a series of weights attributed
% to quality measures of the found ROI clusters. Calls several subfunctions. Saves intermediate mata matrices on disk
% (set of all solutions, matrix and plots of multivariate data, top 5 ranked cluster solutions' topoplot and dipoles).
% If the intermediate steps have already been calculated, you can set the flags to do these steps on 0, then the old
% files will be loaded if available. Stores all relevant information in STUDY.bemobil.
%
% Usage:
%   >>  [STUDY, ALLEEG, EEG] = bemobil_repeated_clustering_and_evaluation(STUDY, ALLEEG, EEG, outlier_sigma, n_clust,...
%       n_iterations, cluster_ROI_talairach, quality_measure_weights, do_clustering, do_multivariate_data,...
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
%   cluster_ROI_talairach           - talairach coordinates of the region of interest (THIS NEEDS TO BE A STRUCT consisting of
%                                   .x, .y, and .z fields)
%   quality_measure_weights         - vector of weights for quality measures. 6 entries: subjects, ICs/subjects, normalized
%                                   spread, mean RV, distance from ROI, mahalanobis distance from median of multivariate
%                                   distribution (put this very high to get the most "normal" solution)
%   do_clustering                   - [0|1] -> whether or not the clustering should be done (it takes a lot of time). If 0, old
%                                   files will be loaded if possible, otherwise will be set to 1 anyways.
%   do_multivariate_data            - [0|1] -> whether or not the creation of the multivariate dataset should be done. If 0, old
%                                   files will be loaded if possible, otherwise will be set to 1 anyways.
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
%   /evaluations/weights_<quality_measure_weights>/rank-<1-5>_ROI_plot.fig topoplots and dipole plots of the top 5
%       solutions
%   
%
% See also:
%   EEGLAB, bemobil_precluster, bemobil_repeated_clustering, bemobil_create_multivariate_data_from_cluster_solutions, bemobil_clustering_rank_solutions
%
% Authors: Marius Klug, 2018


function [STUDY, ALLEEG, EEG] = bemobil_repeated_clustering_and_evaluation(STUDY, ALLEEG, EEG, outlier_sigma, n_clust, n_iterations, cluster_ROI_talairach, quality_measure_weights, do_clustering, do_multivariate_data, filepath_STUDY, filename_STUDY, filepath_clustering_solutions, filename_clustering_solutions, filepath_multivariate_data, filename_multivariate_data)

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

filename_clustering_solutions_with_path = [filepath_clustering_solutions '\' filename_clustering_solutions];

if do_clustering
    
    clustering_solutions = bemobil_repeated_clustering(STUDY,ALLEEG, n_iterations, n_clust, outlier_sigma,STUDY.bemobil.clustering.preclustparams);
    
    disp('Saving clustering solutions...')
    save(filename_clustering_solutions_with_path,'clustering_solutions','-v7.3')
    disp('...done.')
else
    disp('Loading clustering solutions...')
    load(filename_clustering_solutions_with_path,'clustering_solutions')
    disp('...done.')
end

%% create multivariate data of clusters using ROI (takes a little time)

filename_multivariate_data_with_path = [filepath_clustering_solutions '\' filename_multivariate_data];

if do_multivariate_data
    
    [cluster_multivariate_data data_plot] = bemobil_create_multivariate_data_from_cluster_solutions(STUDY,ALLEEG,clustering_solutions,cluster_ROI_talairach);
    
    disp('Saving multivariate data file...')
    save(filename_multivariate_data_with_path,'cluster_multivariate_data','-v7.3')
    savefig(data_plot, [filename_multivariate_data_with_path '_plot'])
    close(data_plot)
    disp('...done.')
    
else
    disp('Loading multivariate data file...')
    load(filename_multivariate_data_with_path,'cluster_multivariate_data')
    disp('...done.')
end

%% find best clustering solutions and store in STUDY. Also crate and save top 5 plots

ranked_solutions = bemobil_clustering_rank_solutions(cluster_multivariate_data,quality_measure_weights);

filepath_to_evaluations = [filepath_clustering_solutions '\evaluations\weights_' num2str(quality_measure_weights) '\'];
mkdir(filepath_to_evaluations)

for rank = [5 4 3 2 1]
    
    STUDY.cluster = clustering_solutions.(['solution_' num2str(ranked_solutions(rank))]);
    
    STUDY.bemobil.clustering.parameters = clustering_solutions.parameters;
    STUDY.bemobil.clustering.cluster_ROI_talairach = cluster_multivariate_data.cluster_ROI_talairach;
    STUDY.bemobil.clustering.cluster_ROI_index = cluster_multivariate_data.best_fitting_cluster(ranked_solutions(rank));
    STUDY.bemobil.clustering.cluster_ROI_multivariate_data = cluster_multivariate_data.data(ranked_solutions(rank),:);
    STUDY.bemobil.clustering.cluster_ROI_multivariate_data
    
    topoplot_and_dipplot_fig = figure;
    subplot(1,2,1);
    std_topoplot(STUDY,ALLEEG,'clusters',STUDY.bemobil.clustering.cluster_ROI_index,'figure','off');
    subplot(1,2,2);
    %     dipplot_fig = figure;
    std_dipplot(STUDY,ALLEEG,'clusters',STUDY.bemobil.clustering.cluster_ROI_index,'figure','off');
    
    savefig(topoplot_and_dipplot_fig, [filepath_to_evaluations 'rank-' num2str(rank) '_ROI_plot'])
end

% find dipole locations, centroids, and residual variances of clusters
STUDY = bemobil_dipoles(STUDY,ALLEEG);
STUDY.bemobil.clustering.quality_measure_weights = quality_measure_weights;
fprintf('Best fitting cluster: %d\n', STUDY.bemobil.clustering.cluster_ROI_index)


% save on disk
disp('Saving STUDY...')
[STUDY, EEG] = pop_savestudy( STUDY, EEG, 'filename',filename_STUDY,'filepath',filepath_STUDY);
disp('...done');
eeglab redraw


%%

%
% for n = 1:5
%     residual_variance_threshold = 1;
% STUDY.cluster = clustering_solutions(ranked_solutions(n)).cluster;
% clustering_solution = STUDY.cluster;
% clusters = best_fitting_cluster(ranked_solutions(n))
%
% cluster_multivariate_data.data(ranked_solutions(n),:)
% standardized_quality_measures(ranked_solutions(n),:)
% std_dipplot(STUDY, ALLEEG, 'clusters', clusters);
% plotERSPs_2
% end
% disp('...done')
% %% plot for evaluating if selected cluster is the right one
%
% % load(['P:\Lukas_Gehrke\studies\Spot_Rotation\data\5_study_level\clustering_solutions\clustering-solutions_' num2str(n_clust) '-cluster_' num2str(sigma) '-sigma_' num2str(dipole_weight) '-dipoles_' num2str(spec_weight) '-spectra_' num2str(ersp_weight) '-ersp_' num2str(max_iterations) '-iterations'],'clustering_solutions')
% % solutions = 1:length(fields(clustering_solutions))';
% solutions = closest_median_solutions';
%
% for solution = solutions
%     solution
%     STUDY.cluster = clustering_solutions(solution).cluster;
% %
% %     num2str(best_fitting_cluster_distance(solution))
% %     num2str(best_fitting_cluster_n_subjects(solution))
% %     num2str(best_fitting_cluster_n_ICs(solution))
% %     num2str(best_fitting_cluster_spread(solution))
% %     num2str(best_fitting_cluster_median_rv(solution))
% %     num2str(best_fitting_cluster_mean_rv(solution))
% %     num2str(best_fitting_cluster_normalized_spread(solution))
% % best_fitting_cluster_distance(solution)
% % best_fitting_cluster_n_subjects(solution)
% % best_fitting_cluster_n_ICs(solution)
% % best_fitting_cluster_spread(solution)
% % best_fitting_cluster_median_rv(solution)
% % best_fitting_cluster_mean_rv(solution)
% % best_fitting_cluster_normalized_spread(solution)
%
% %     f = figure;
%     std_dipplot(STUDY,ALLEEG,'clusters',best_fitting_cluster(solution));%,'figure','off');
%     g = figure;
%     title(['BEST FITTING CLUSTER: ' num2str(best_fitting_cluster(solution))...
%     ', distance: ' num2str(best_fitting_cluster_distance(solution))...
%     ', mean rv: ' num2str(best_fitting_cluster_mean_rv(solution))...
%     ', normalized spread: ' num2str(best_fitting_cluster_normalized_spread(solution))...
%     ]);
%     std_topoplot(STUDY,ALLEEG,'clusters',3:length(STUDY.cluster),'figure','off');
%
% %     button = waitforbuttonpress;
% %     if button == 1
% %         close(f)
% %         close(g)
% %     end
% end