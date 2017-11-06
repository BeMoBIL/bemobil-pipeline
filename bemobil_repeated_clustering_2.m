freqrange = [3 30];
dipole_weight = 10;
ersp_weight = 0;
spec_weight = 1;
sigma = 3;
n_clust = 50;
max_iterations = 200;
iterations = 1:max_iterations;


% [STUDY ALLEEG] = std_preclust(STUDY, ALLEEG, 1,{'spec' 'npca' 10 'norm' 1 'weight' spec_weight 'freqrange' freqrange },{'dipoles' 'norm' 1 'weight' dipole_weight},{'ersp' 'npca' 10 'freqrange' freqrange 'timewindow' [] 'norm' 1 'weight' ersp_weight},{'finaldim' 'npca' 10});
[STUDY ALLEEG] = std_preclust(STUDY, ALLEEG, 1,{'spec' 'npca' 10 'norm' 1 'weight' spec_weight 'freqrange' freqrange },{'dipoles' 'norm' 1 'weight' dipole_weight},{'finaldim' 'npca' 10});


%% cluster

for iteration = iterations
    
    iteration
    % clustering
    [STUDY] = pop_clust(STUDY, ALLEEG, 'algorithm','kmeans','clus_num',  n_clust , 'outliers',  sigma );
    
    % create the average topoplots for the clusters. cluster 1 is the
    % parentcluster, 2 are the outliers. also store dipole (and centroid) locations
    [STUDY, centroid] = std_readtopoclust(STUDY,ALLEEG, 3:length(STUDY.cluster));
    STUDY = bemobil_dipoles(STUDY,ALLEEG);
    
    % save the clustering results
    clustering_solutions.(['solution_' num2str(iteration)]) = STUDY.cluster;
    
end

save(['P:\Lukas_Gehrke\studies\Spot_Rotation\data\5_study_level\clustering_solutions\clustering-solutions_' num2str(n_clust) '-cluster_' num2str(sigma) '-sigma_' num2str(dipole_weight) '-dipoles_' num2str(spec_weight) '-spectra_' num2str(ersp_weight) '-ersp_' num2str(max_iterations) '-iterations'],'clustering_solutions','-v7.3')

%% plot 

for iteration = iterations

STUDY.cluster = clustering_solutions.(['solution_' num2str(iteration)]);
std_topoplot(STUDY,ALLEEG,'clusters',3:length(STUDY.cluster));

end