freqrange = [4 30];
dipole_weight = 6;
ersp_weight = 3;
ersp_timewindow = [0 6326];
spec_weight = 0;
scalp_weight = 1;

sigma = 3;
n_clust = 50;
max_iterations = 10000;
iterations = 1:max_iterations;
cluster_ROI_talairach = [0 -45 10];

%% precluster

% [STUDY ALLEEG] = std_preclust(STUDY, ALLEEG, 1,{'spec' 'npca' 10 'norm' 1 'weight' spec_weight 'freqrange' freqrange },{'dipoles' 'norm' 1 'weight' dipole_weight},{'ersp' 'npca' 10 'freqrange' freqrange 'timewindow' ersp_timewindow 'norm' 1 'weight' ersp_weight},{'finaldim' 'npca' 10});
% [STUDY ALLEEG] = std_preclust(STUDY, ALLEEG, 1,{'spec' 'npca' 10 'norm' 1 'weight' spec_weight 'freqrange' freqrange },{'dipoles' 'norm' 1 'weight' dipole_weight},{'finaldim' 'npca' 10});
% [STUDY ALLEEG] = std_preclust(STUDY, ALLEEG, 1,{'dipoles' 'norm' 1 'weight' dipole_weight},{'ersp' 'npca' 10 'freqrange' freqrange 'timewindow' ersp_timewindow 'norm' 1 'weight' ersp_weight},{'finaldim' 'npca' 10});
[STUDY ALLEEG] = std_preclust(STUDY, ALLEEG, 1,{'scalp' 'npca' 10 'norm' 1 'weight' scalp_weight 'abso' 1},{'dipoles' 'norm' 1 'weight' dipole_weight},{'ersp' 'npca' 10 'freqrange' freqrange 'timewindow' ersp_timewindow  'norm' 1 'weight' ersp_weight},{'finaldim' 'npca' 10});

%% cluster

for iteration = iterations
    
    disp(['Iteration: ' num2str(iteration)])
    % clustering
    [STUDY] = pop_clust(STUDY, ALLEEG, 'algorithm','kmeans','clus_num',  n_clust , 'outliers',  sigma );
    
    % create the average topoplots for the clusters. cluster 1 is the
    % parentcluster, 2 are the outliers. 
%     [STUDY, centroid] = std_readtopoclust(STUDY,ALLEEG, 3:length(STUDY.cluster));
    
    % save the clustering results
    clustering_solutions.(['solution_' num2str(iteration)]) = STUDY.cluster;
    
end
disp('Saving clustering solutions...')
save(['P:\Lukas_Gehrke\studies\Spot_Rotation\data\5_study_level\clustering_solutions\clustering-solutions_' num2str(n_clust) '-cluster_' num2str(sigma) '-sigma_' num2str(dipole_weight) '-dipoles_' num2str(scalp_weight) '-scalp_' num2str(ersp_weight) '-ersp_' num2str(max_iterations) '-iterations'],'clustering_solutions','-v7.3')
disp('...done!')
%% evaluate clusters

disp('Loading clustering solutions...')
% load(['P:\Lukas_Gehrke\studies\Spot_Rotation\data\5_study_level\clustering_solutions\clustering-solutions_' num2str(n_clust) '-cluster_' num2str(sigma) '-sigma_' num2str(dipole_weight) '-dipoles_' num2str(scalp_weight) '-scalp_' num2str(ersp_weight) '-ersp_' num2str(max_iterations) '-iterations'],'clustering_solutions')
disp('...done')
best_fitting_cluster = zeros(length(fields(clustering_solutions)),1);
best_fitting_cluster_distance = inf(length(fields(clustering_solutions)),1);
best_fitting_cluster_x = inf(length(fields(clustering_solutions)),1);
best_fitting_cluster_y = inf(length(fields(clustering_solutions)),1);
best_fitting_cluster_z = inf(length(fields(clustering_solutions)),1);
best_fitting_cluster_spread = inf(length(fields(clustering_solutions)),1);
best_fitting_cluster_n_subjects = inf(length(fields(clustering_solutions)),1);
best_fitting_cluster_n_ICs = inf(length(fields(clustering_solutions)),1);
best_fitting_cluster_normalized_spread = inf(length(fields(clustering_solutions)),1);
best_fitting_cluster_median_rv = inf(length(fields(clustering_solutions)),1);
best_fitting_cluster_mean_rv = inf(length(fields(clustering_solutions)),1);
d_median = inf(length(fields(clustering_solutions)),1);

for solution = 1:length(fields(clustering_solutions))
    solution
    
    STUDY.cluster = clustering_solutions.(['solution_' num2str(solution)]);
    
    % find dipole locations and centroids
    STUDY = bemobil_dipoles(STUDY,ALLEEG);

    % find closest cluster to ROI
    for cluster = 3:length(STUDY.cluster)

        tmp_distance = sqrt((STUDY.cluster(cluster).dipole.posxyz(1) - cluster_ROI_talairach(1))^2 + (STUDY.cluster(cluster).dipole.posxyz(2) - cluster_ROI_talairach(2))^2 + (STUDY.cluster(cluster).dipole.posxyz(3) - cluster_ROI_talairach(3))^2);

        if tmp_distance < best_fitting_cluster_distance(solution)
            % declare this cluster to be the best fitting cluster of this solution
            best_fitting_cluster(solution) = cluster;
            
            % take quality measures
            best_fitting_cluster_distance(solution) = tmp_distance;
            best_fitting_cluster_x(solution) = STUDY.cluster(cluster).dipole.posxyz(1);
            best_fitting_cluster_y(solution) = STUDY.cluster(cluster).dipole.posxyz(2);
            best_fitting_cluster_z(solution) = STUDY.cluster(cluster).dipole.posxyz(3);
            best_fitting_cluster_n_subjects(solution) = length(unique(STUDY.cluster(cluster).sets));
            best_fitting_cluster_n_ICs(solution) = length(STUDY.cluster(cluster).comps);
            best_fitting_cluster_spread(solution) = STUDY.cluster(cluster).squared_deviations;
            best_fitting_cluster_median_rv(solution) = STUDY.cluster(cluster).median_rv;
            best_fitting_cluster_mean_rv(solution) = STUDY.cluster(cluster).mean_rv;
            best_fitting_cluster_normalized_spread(solution) = best_fitting_cluster_spread(solution) / best_fitting_cluster_n_ICs(solution);
            
        end

    end

end

%%
clear standardized_quality_measures

% plot quality measures
% figure;hist(best_fitting_cluster_n_subjects); title(['best fitting cluster number of subjects, median = ' num2str(median(best_fitting_cluster_n_subjects))])
% figure;hist(best_fitting_cluster_n_ICs); title(['best fitting cluster number of ICs, median = ' num2str(median(best_fitting_cluster_n_ICs))])
% figure;hist(best_fitting_cluster_spread); title(['best fitting cluster spread, median = ' num2str(median(best_fitting_cluster_spread))])
% figure;hist(best_fitting_cluster_normalized_spread); title(['best fitting cluster normalized spread, median = ' num2str(median(best_fitting_cluster_normalized_spread))])
% figure;hist(best_fitting_cluster_distance); title(['best fitting cluster distance, median = ' num2str(median(best_fitting_cluster_distance))])
% figure;hist(best_fitting_cluster_x); title(['best fitting cluster x, median = ' num2str(median(best_fitting_cluster_x))])
% figure;hist(best_fitting_cluster_y); title(['best fitting cluster y, median = ' num2str(median(best_fitting_cluster_y))])
% figure;hist(best_fitting_cluster_z); title(['best fitting cluster z, median = ' num2str(median(best_fitting_cluster_z))])
% figure;hist(best_fitting_cluster_median_rv); title(['best fitting cluster median residual variance, median = ' num2str(median(best_fitting_cluster_median_rv))])
% figure;hist(best_fitting_cluster_mean_rv); title(['best fitting cluster mean residual variance, median = ' num2str(median(best_fitting_cluster_mean_rv))])

% create multivariate statistics data matrix -> non-normalized spread will not be used
clustering_solutions_multivariate_data = zeros(length(fields(clustering_solutions)),4);
% clustering_solutions_multivariate_data(:,1) = best_fitting_cluster_distance;
clustering_solutions_multivariate_data(:,1) = best_fitting_cluster_n_subjects;
clustering_solutions_multivariate_data(:,2) = best_fitting_cluster_n_ICs;
clustering_solutions_multivariate_data(:,3) = best_fitting_cluster_normalized_spread;
clustering_solutions_multivariate_data(:,4) = best_fitting_cluster_mean_rv;
clustering_solutions_multivariate_data(:,5) = best_fitting_cluster_x;
clustering_solutions_multivariate_data(:,6) = best_fitting_cluster_y;
clustering_solutions_multivariate_data(:,7) = best_fitting_cluster_z;

clustering_solutions_multivariate_median = median(clustering_solutions_multivariate_data,1);
clustering_solutions_multivariate_covariance = cov(clustering_solutions_multivariate_data);

% find mahalanobis distance from median for each solution
d_mean = mahal(clustering_solutions_multivariate_data,clustering_solutions_multivariate_data);
for solution = 1:length(fields(clustering_solutions))
   
    d_median(solution) = (clustering_solutions_multivariate_data(solution,:)-clustering_solutions_multivariate_median)*inv(clustering_solutions_multivariate_covariance)*(clustering_solutions_multivariate_data(solution,:)-clustering_solutions_multivariate_median)';
    
end

n = 10;
[d_median_s, index] = sort(d_median);
closest_median_solutions = index(1:n);

[d_mean_s, index] = sort(d_mean);
closest_mean_solutions = index(1:n);


n=1;
STUDY.cluster = clustering_solutions.(['solution_' num2str(closest_median_solutions(n))]);
clustering_solution = STUDY.cluster;
closest_median_solutions(n);
clusters = best_fitting_cluster(closest_median_solutions(n));
clustering_solutions_multivariate_data(closest_median_solutions(n),:);

% alternative method: weigh standardized factors and choose the best solution
standardized_quality_measures(:,1) = best_fitting_cluster_n_subjects ./ max(best_fitting_cluster_n_subjects); % should be high
standardized_quality_measures(:,2) = best_fitting_cluster_distance ./ max(best_fitting_cluster_distance); % should be low
standardized_quality_measures(:,3) = d_median ./ max(d_median); % should be low
standardized_quality_measures(:,4) = (best_fitting_cluster_n_ICs./best_fitting_cluster_n_subjects)./ max(best_fitting_cluster_n_ICs./best_fitting_cluster_n_subjects); % should be low
% standardized_quality_measures(:,4) = best_fitting_cluster_n_ICs ./ max(best_fitting_cluster_n_ICs); % should be low
standardized_quality_measures(:,5) = best_fitting_cluster_normalized_spread ./ max(best_fitting_cluster_normalized_spread); % should be low
standardized_quality_measures(:,6) = best_fitting_cluster_mean_rv ./ max(best_fitting_cluster_mean_rv); % should be low


weights=[2,-3,-1,-1,-1,-2]; % subjects, dist from ROI, mahal dist from median, #ICs/#S ratio, spread, rv
weighted_scores = standardized_quality_measures .* weights;
sum_of_weighted_scores = sum(weighted_scores,2);

[sum_of_weighted_scores_sorted, best_solutions] = sort(sum_of_weighted_scores,'descend');
%%
for n = 1:5
    residual_variance_threshold = 1;
STUDY.cluster = clustering_solutions.(['solution_' num2str(best_solutions(n))]);
clustering_solution = STUDY.cluster;
clusters = best_fitting_cluster(best_solutions(n))

clustering_solutions_multivariate_data(best_solutions(n),:)
standardized_quality_measures(best_solutions(n),:)
std_dipplot(STUDY, ALLEEG, 'clusters', clusters);
plotERSPs_2
end
disp('...done')
%% plot for evaluating if selected cluster is the right one

% load(['P:\Lukas_Gehrke\studies\Spot_Rotation\data\5_study_level\clustering_solutions\clustering-solutions_' num2str(n_clust) '-cluster_' num2str(sigma) '-sigma_' num2str(dipole_weight) '-dipoles_' num2str(spec_weight) '-spectra_' num2str(ersp_weight) '-ersp_' num2str(max_iterations) '-iterations'],'clustering_solutions')
% solutions = 1:length(fields(clustering_solutions))';
solutions = closest_median_solutions';

for solution = solutions
    solution
    STUDY.cluster = clustering_solutions.(['solution_' num2str(solution)]);
%     
%     num2str(best_fitting_cluster_distance(solution))
%     num2str(best_fitting_cluster_n_subjects(solution))
%     num2str(best_fitting_cluster_n_ICs(solution))
%     num2str(best_fitting_cluster_spread(solution))
%     num2str(best_fitting_cluster_median_rv(solution))
%     num2str(best_fitting_cluster_mean_rv(solution))
%     num2str(best_fitting_cluster_normalized_spread(solution))
% best_fitting_cluster_distance(solution)
% best_fitting_cluster_n_subjects(solution)
% best_fitting_cluster_n_ICs(solution)
% best_fitting_cluster_spread(solution)
% best_fitting_cluster_median_rv(solution)
% best_fitting_cluster_mean_rv(solution)
% best_fitting_cluster_normalized_spread(solution)
    
%     f = figure;
    std_dipplot(STUDY,ALLEEG,'clusters',best_fitting_cluster(solution));%,'figure','off');
    g = figure;
    title(['BEST FITTING CLUSTER: ' num2str(best_fitting_cluster(solution))...
    ', distance: ' num2str(best_fitting_cluster_distance(solution))...
    ', mean rv: ' num2str(best_fitting_cluster_mean_rv(solution))...
    ', normalized spread: ' num2str(best_fitting_cluster_normalized_spread(solution))...
    ]);
    std_topoplot(STUDY,ALLEEG,'clusters',3:length(STUDY.cluster),'figure','off');
    
%     button = waitforbuttonpress;
%     if button == 1
%         close(f)
%         close(g)
%     end
end