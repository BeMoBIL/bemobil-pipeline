function clustering_solutions = bemobil_repeated_clustering_full(STUDY,ALLEEG, n_iterations, n_clust, outlier_sigma,preclustparams)

fprintf('Clustering %d times...\n', n_iterations)
    
    clustering_solutions = struct;
    
    for iteration = 1:n_iterations
        tic
        
        fprintf('Iteration %d of %d \n', iteration,n_iterations)
        % clustering
        [STUDY] = pop_clust(STUDY, ALLEEG, 'algorithm','kmeans','clus_num',  n_clust , 'outliers',  outlier_sigma );
        
        % create the average topoplots for the clusters. cluster 1 is the parentcluster, 2 are the
        % outliers. For testing purposes only.
        % [STUDY, centroid] = std_readtopoclust(STUDY,ALLEEG, 3:length(STUDY.cluster));
        
        clustering_solutions.(['solution_' num2str(iteration)]) = STUDY.cluster;
        
        lastduration = toc; % stop checking the time and plot estimated time of arrival
        eta = lastduration * (n_iterations-iteration);
        fprintf('Last duration: %1.2f s \n',round(lastduration,2))
        fprintf('ETA: %d h, %d min \n', floor(eta/3600), round(((eta/3600)-floor(eta/3600))*60))
        
    end
    
    % store clustering info
    clustering_solutions.parameters.preclustparams = preclustparams;
    clustering_solutions.parameters.outlier_sigma = outlier_sigma;
    clustering_solutions.parameters.n_clust = n_clust;
    clustering_solutions.parameters.n_iterations = n_iterations;