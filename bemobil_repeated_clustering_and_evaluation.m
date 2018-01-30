function [STUDY, ALLEEG, EEG] = bemobil_repeated_clustering_and_evaluation(STUDY, ALLEEG, EEG, outlier_sigma, n_clust, n_iterations, cluster_ROI_talairach, quality_measure_weights, do_clustering, do_multivariate_data, filepath_STUDY, filename_STUDY, filepath_clustering_solutions, filename_clustering_solutions, filepath_multivariate_data, filename_multivariate_data)

% only save a STUDY file on disk if both a name and a path are provided
save_file_on_disk = exist('filename_STUDY', 'var');

% check if files already exist and show warning if it does
% STUDY
if save_file_on_disk
    mkdir(filepath_STUDY); % make sure that folder exists, nothing happens if so
    dir_files = dir(filepath_STUDY);
    if ismember(filename_STUDY, {dir_files.name})
        warning([filename_STUDY ' file already exists in: ' filepath_STUDY '. File will be overwritten...']);
    end
end

% clustering_solutions
mkdir(filepath_clustering_solutions); % make sure that folder exists, nothing happens if so
dir_files = dir(filepath_clustering_solutions);

if ismember(filename_clustering_solutions, {dir_files.name})
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

if ismember(filename_multivariate_data, {dir_files.name})
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

filename_multivariate_data_with_path = [filepath_clustering_solutions '\' filename_clustering_solutions];

if do_multivariate_data

    [cluster_multivariate_data data_plot] = bemobil_create_multivariate_data_from_cluster_solutions(STUDY,ALLEEG,clustering_solutions,cluster_ROI_talairach);
    
    disp('Saving multivariate data file...')
    save(filename_multivariate_data_with_path,'cluster_multivariate_data','-v7.3')
    savefig(data_plot, [filename_multivariate_data_with_path '_' data_plot])
    disp('...done.')
    
else
    disp('Loading multivariate data file...')
    load(filename_clustering_solutions_with_path,'cluster_multivariate_data')
    disp('...done.')
end

%% find best clustering solutions and store in STUDY

ranked_solutions = bemobil_clustering_rank_solutions(cluster_multivariate_data,quality_measure_weights);

STUDY.cluster = clustering_solutions(ranked_solutions(1)).cluster;

STUDY.bemobil.clustering.parameters = clustering_solutions.parameters;
STUDY.bemobil.clustering.cluster_ROI_talairach = cluster_multivariate_data.cluster_ROI_talairach;
STUDY.bemobil.clustering.cluster_ROI_index = cluster_multivariate_data.best_fitting_cluster(ranked_solutions(1));

% save on disk
if save_file_on_disk
    [STUDY, EEG] = pop_savestudy( STUDY, EEG, 'filename',STUDY_filename,'filepath',STUDY_filepath);
    disp('...done');
end

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