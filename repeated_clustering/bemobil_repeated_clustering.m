% bemobil_repeated_clustering() - Repeatedly clusters the dipoles of the EEGLAB STUDY data set according to parameters.
% Saves relevant parameter information in clustering_solutions.parameters. The preclustering should have been done using
% bemobil_precluster()
% 
%
% Usage:
%   >>  clustering_solutions = bemobil_repeated_clustering(STUDY,ALLEEG, n_iterations, n_clust, outlier_sigma,preclustparams)
%   
% Inputs:
%   ALLEEG                          - complete EEGLAB data set structure
%   STUDY                           - STUDY data set of EEGLAB, which has to be loaded previously
%   outlier_sigma                   - standard deviations boundary for outlier detection (e.g. 3)
%   n_clust                         - number of clusters to be created 
%   n_iterations                    - number of iterations to be performed for repeated clustering
%   preclustparams                  - parameters which have been used for preclustering. Take them from 
%                                   STUDY.bemobil.clustering.preclustparams
%   
%
% Outputs:
%   clustering_solutions            - struct of all clustering solutions to be evaluated later
%  
%
% See also:
%   EEGLAB, bemobil_precluster, bemobil_repeated_clustering_and_evaluation
%
% Authors: Marius Klug, 2018

function clustering_solutions = bemobil_repeated_clustering(STUDY,ALLEEG, n_iterations, n_clust, outlier_sigma,preclustparams)

fprintf('Clustering %d times...\n', n_iterations)
    
    clustering_solutions = struct;
    
    for iteration = 1:n_iterations
        
        % start time
        tic
        
        fprintf('Iteration %d of %d \n', iteration,n_iterations)
        
        % clustering according to parameters
        [STUDY] = pop_clust(STUDY, ALLEEG, 'algorithm','kmeans','clus_num',  n_clust , 'outliers',  outlier_sigma );
        
        % create the average topoplots for the clusters. cluster 1 is the parentcluster, 2 are the
        % outliers. For testing purposes only.
        % [STUDY, centroid] = std_readtopoclust(STUDY,ALLEEG, 3:length(STUDY.cluster));
        
        % store info in output struct
        clustering_solutions.(['solution_' num2str(iteration)]) = STUDY.cluster;
        
        % stop checking the time and plot estimated time of arrival
        lastduration = toc; 
        eta = lastduration * (n_iterations-iteration);
        fprintf('Last duration: %1.2f s \n',round(lastduration,2))
        fprintf('ETA: %d h, %d min \n', floor(eta/3600), round(((eta/3600)-floor(eta/3600))*60))
        
    end
    
    % store clustering info
    clustering_solutions.parameters.preclustparams = preclustparams;
    clustering_solutions.parameters.outlier_sigma = outlier_sigma;
    clustering_solutions.parameters.n_clust = n_clust;
    clustering_solutions.parameters.n_iterations = n_iterations;