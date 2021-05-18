% bemobil_clustering_rank_solutions() - Ranks the complete multivariate data set of repeated clusterings according to a
% given set of weights for the quality measures derived from the multivariate data
% 
% Usage:
%   >>  ranked_solutions = bemobil_clustering_rank_solutions(cluster_multivariate_data,quality_measure_weights)
%
% Inputs:
%   cluster_multivariate_data   - multivariate data set for each solutions' best fitting cluster (to ROI). Struct with
%                               fields data (the data set consisting of these values: n_subjects,n_ICs,n_ICs/n_subjects,...
%                               normalized_spread,mean_rv,x,y,z,distance_from_ROI), best_fitting_cluster (the number of
%                               the best fitting cluster in the solution), and cluster_ROI_MNI (the input)
%   quality_measure_weights     - vector of weights for quality measures. 6 entries: subjects, ICs/subjects, normalized
%                               spread, mean RV, distance from ROI, mahalanobis distance from median of multivariate
%                               distribution (put this very high to get the most "normal" solution)
%
% Outputs:
%   
%   ranked_solutions            - the solutions sorted according to their rank. ranked_solutions(1) is the best
%                               solution's number in the clustering_solutions data set according to the given weights.
%   
%
% See also:
%   EEGLAB, bemobil_repeated_clustering_and_evaluation, bemobil_repeated_clustering, bemobil_create_multivariate_data_from_cluster_solutions
%
% Authors: Marius Klug, 2018

function ranked_solutions = bemobil_clustering_rank_solutions(cluster_multivariate_data,quality_measure_weights)

clustering_solutions_multivariate_median = median(cluster_multivariate_data.data,1);
clustering_solutions_multivariate_covariance = cov(cluster_multivariate_data.data);

% find mahalanobis distance from median for each solution
% d_mean = mahal(cluster_multivariate_data.data,cluster_multivariate_data.data);
for solution = 1:size(cluster_multivariate_data.data,1)
   
    d_median(solution) = (cluster_multivariate_data.data(solution,:)-clustering_solutions_multivariate_median)*inv(clustering_solutions_multivariate_covariance)*(cluster_multivariate_data.data(solution,:)-clustering_solutions_multivariate_median)';
    
end

% create standardized measures
standardized_quality_measures(:,1) = cluster_multivariate_data.data(:,1) ./ max(cluster_multivariate_data.data(:,1)); % subjects
standardized_quality_measures(:,2) = cluster_multivariate_data.data(:,3) ./ max(cluster_multivariate_data.data(:,3)); % ICs/subjects
standardized_quality_measures(:,3) = cluster_multivariate_data.data(:,4) ./ max(cluster_multivariate_data.data(:,4)); % normalized spread
standardized_quality_measures(:,4) = cluster_multivariate_data.data(:,5) ./ max(cluster_multivariate_data.data(:,5)); % mean RV
standardized_quality_measures(:,5) = cluster_multivariate_data.data(:,9) ./ max(cluster_multivariate_data.data(:,9)); % distance from ROI
standardized_quality_measures(:,6) = d_median ./ max(d_median);                                                       % mahal distance from median

weighted_scores = standardized_quality_measures .* quality_measure_weights;
sum_of_weighted_scores = sum(weighted_scores,2);
[~, ranked_solutions] = sort(sum_of_weighted_scores,'descend');