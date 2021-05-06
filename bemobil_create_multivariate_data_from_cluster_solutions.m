% bemobil_create_multivariate_data_from_cluster_solutions() - Creates a multivariate data set with all available info of
% each cluster solution to be used for evaluation later. A region of interest is necessary: Each clustering solution's
% closes cluster to that ROI is being evaluated by finding it's specific values.
% 
% Usage:
%   >>  [cluster_multivariate_data, data_plot] = bemobil_create_multivariate_data_from_cluster_solutions(STUDY,ALLEEG,...
%       clustering_solutions,cluster_ROI_MNI)
%
% Inputs:
%   ALLEEG                      - complete EEGLAB data set structure
%   STUDY                       - STUDY data set of EEGLAB, which has to be loaded previously
%   clustering_solutions        - struct of all clustering solutions created from bemobil_repeated_clustering
%   cluster_ROI_MNI             - MNI coordinates of the region of interest (THIS NEEDS TO BE A STRUCT consisting of
%                                   .x, .y, and .z fields)
%
% Outputs:
%   cluster_multivariate_data   - multivariate data set for each solutions' best fitting cluster (to ROI). Struct with
%                               fields data (the data set consisting of these values: n_subjects,n_ICs,n_ICs/n_subjects,...
%                               normalized_spread,mean_rv,x,y,z,distance_from_ROI), best_fitting_cluster (the number of
%                               the best fitting cluster in the solution), and cluster_ROI_MNI (the input)
%   data_plot                   - plot of the histograms of all dimensions
%   
%
% See also:
%   EEGLAB, bemobil_repeated_clustering_and_evaluation, bemobil_repeated_clustering
%
% Authors: Marius Klug, 2018

function [cluster_multivariate_data, data_plot] = bemobil_create_multivariate_data_from_cluster_solutions(STUDY,ALLEEG,clustering_solutions,cluster_ROI_MNI)

best_fitting_cluster = zeros((length(fields(clustering_solutions))-1),1);
best_fitting_cluster_distance = inf((length(fields(clustering_solutions))-1),1);
best_fitting_cluster_x = inf((length(fields(clustering_solutions))-1),1);
best_fitting_cluster_y = inf((length(fields(clustering_solutions))-1),1);
best_fitting_cluster_z = inf((length(fields(clustering_solutions))-1),1);
best_fitting_cluster_spread = inf((length(fields(clustering_solutions))-1),1);
best_fitting_cluster_n_subjects = inf((length(fields(clustering_solutions))-1),1);
best_fitting_cluster_n_ICs = inf((length(fields(clustering_solutions))-1),1);
best_fitting_cluster_normalized_spread = inf((length(fields(clustering_solutions))-1),1);
best_fitting_cluster_mean_rv = inf((length(fields(clustering_solutions))-1),1);
d_median = inf((length(fields(clustering_solutions))-1),1);

tic % start time measure
for solution = 1:(length(fields(clustering_solutions))-1)
    
    STUDY.cluster = clustering_solutions.(['solution_' num2str(solution)]);
    
    % find dipole locations and centroids
    STUDY = bemobil_dipoles(STUDY,ALLEEG);
    
    % find closest cluster to ROI
    for cluster = 3:length(STUDY.cluster) % cluster 1 contains all ICs and cluster 2 contains the outliers
        
		% Pythagoras in 3D
        tmp_distance = sqrt((STUDY.cluster(cluster).dipole.posxyz(1) - cluster_ROI_MNI.x)^2 +...
							(STUDY.cluster(cluster).dipole.posxyz(2) - cluster_ROI_MNI.y)^2 +...
							(STUDY.cluster(cluster).dipole.posxyz(3) - cluster_ROI_MNI.z)^2);
        
        if tmp_distance < best_fitting_cluster_distance(solution)
            % declare this cluster to be the best fitting cluster of this solution
            best_fitting_cluster(solution) = cluster;
            
            % take quality measures
            best_fitting_cluster_distance(solution) = tmp_distance;
            best_fitting_cluster_x(solution) = STUDY.cluster(cluster).dipole.posxyz(1);
            best_fitting_cluster_y(solution) = STUDY.cluster(cluster).dipole.posxyz(2);
            best_fitting_cluster_z(solution) = STUDY.cluster(cluster).dipole.posxyz(3);
            best_fitting_cluster_n_subjects(solution) = length(unique(STUDY.cluster(cluster).sets(1,:)));
            best_fitting_cluster_n_ICs(solution) = length(STUDY.cluster(cluster).comps);
            best_fitting_cluster_spread(solution) = STUDY.cluster(cluster).squared_deviations;
            best_fitting_cluster_mean_rv(solution) = STUDY.cluster(cluster).mean_rv;
            best_fitting_cluster_normalized_spread(solution) = best_fitting_cluster_spread(solution) / best_fitting_cluster_n_ICs(solution);
            
        end
        
    end
    if mod(solution,100) == 0 % print current status every 100 solutions
        fprintf('Solution %d of %d \n', solution,(length(fields(clustering_solutions))-1))
        lastduration = toc; % stop checking the time and plot estimated time of arrival
        eta = lastduration * ((length(fields(clustering_solutions))-1)-solution)/100;
        fprintf('Last duration: %3.2f s \n',round(lastduration,2))
        fprintf('ETA: %d h, %d min \n', floor(eta/3600), round(((eta/3600)-floor(eta/3600))*60))
        tic
    end
    
end

% best_fitting_clusters = struct([]);
% best_fitting_clusters.distance = best_fitting_cluster_distance;
% best_fitting_clusters.x = best_fitting_cluster_x;
% best_fitting_clusters.y = best_fitting_cluster_y;
% best_fitting_clusters.z = best_fitting_cluster_z;
% best_fitting_clusters.n_subjects = best_fitting_cluster_n_subjects;
% best_fitting_clusters.n_ICs = best_fitting_cluster_n_ICs;
% best_fitting_clusters.best_fitting_cluster_spread = best_fitting_cluster_best_fitting_cluster_spread;
% best_fitting_clusters.normalized_spread = best_fitting_cluster_normalized_spread;
% best_fitting_clusters.mean_rv = best_fitting_cluster_mean_rv;

% plot data
data_plot = figure;
subplot(3,3,1);histogram(best_fitting_cluster_n_subjects); title(['best fitting cluster number of subjects, median = ' num2str(median(best_fitting_cluster_n_subjects))])
subplot(3,3,2);histogram(best_fitting_cluster_n_ICs); title(['best fitting cluster number of ICs, median = ' num2str(median(best_fitting_cluster_n_ICs))])
subplot(3,3,3);histogram(best_fitting_cluster_spread); title(['best fitting cluster spread, median = ' num2str(median(best_fitting_cluster_spread))])
subplot(3,3,4);histogram(best_fitting_cluster_normalized_spread); title(['best fitting cluster normalized spread, median = ' num2str(median(best_fitting_cluster_normalized_spread))])
subplot(3,3,5);histogram(best_fitting_cluster_distance); title(['best fitting cluster distance, median = ' num2str(median(best_fitting_cluster_distance))])
subplot(3,3,6);histogram(best_fitting_cluster_mean_rv); title(['best fitting cluster mean residual variance, median = ' num2str(median(best_fitting_cluster_mean_rv))])
subplot(3,3,7);histogram(best_fitting_cluster_x); title(['best fitting cluster x, median = ' num2str(median(best_fitting_cluster_x))])
subplot(3,3,8);histogram(best_fitting_cluster_y); title(['best fitting cluster y, median = ' num2str(median(best_fitting_cluster_y))])
subplot(3,3,9);histogram(best_fitting_cluster_z); title(['best fitting cluster z, median = ' num2str(median(best_fitting_cluster_z))])

% create multivariate statistics data matrix -> non-normalized spread will not be used because it's
% biased towards fewer ICs
multivariate_data = zeros((length(fields(clustering_solutions))-1),4);



multivariate_data(:,1) = best_fitting_cluster_n_subjects;
multivariate_data(:,2) = best_fitting_cluster_n_ICs;
multivariate_data(:,3) = best_fitting_cluster_n_ICs./best_fitting_cluster_n_subjects;
multivariate_data(:,4) = best_fitting_cluster_normalized_spread;
multivariate_data(:,5) = best_fitting_cluster_mean_rv;
multivariate_data(:,6) = best_fitting_cluster_x;
multivariate_data(:,7) = best_fitting_cluster_y;
multivariate_data(:,8) = best_fitting_cluster_z;
multivariate_data(:,9) = best_fitting_cluster_distance;


cluster_multivariate_data.data = multivariate_data;
cluster_multivariate_data.best_fitting_cluster = best_fitting_cluster;
cluster_multivariate_data.cluster_ROI_MNI = cluster_ROI_MNI;