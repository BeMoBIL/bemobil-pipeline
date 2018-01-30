%% load and prepare the stuff
% 
% [ALLEEG EEG CURRENTSET ALLCOM] = eeglab;
% 
% filename = 'outward_rotation_short.study';
% filepath = 'P:\Lukas_Gehrke\studies\Spot_Rotation\data\5_study_level';
% 
% [STUDY ALLEEG] = pop_loadstudy('filename', filename, 'filepath', filepath);
% CURRENTSTUDY = 1; EEG = ALLEEG; CURRENTSET = [1:length(EEG)];
% eeglab('redraw');
% 
% [STUDY ALLEEG] = std_preclust(STUDY, ALLEEG, 1,{'dipoles' 'norm' 1 'weight' 2},{'ersp' 'npca' 10 'freqrange' [] 'timewindow' [] 'norm' 1 'weight' 1});


%% do the stuff
% all_clusters = struct();
max_iterations = 100;
iterations = 1:max_iterations;

for iteration = iterations
    
    iteration
    % clustering
    [STUDY] = pop_clust(STUDY, ALLEEG, 'algorithm','kmeans','clus_num',  80 , 'outliers',  3 );
    
    % create the average topoplots for the clusters. cluster 1 is the
    % parentcluster, 2 are the outliers
    [STUDY, centroid] = std_readtopoclust(STUDY,ALLEEG, 3:length(STUDY.cluster));
    
    % replace NaN values with 0, for correlation later
    for clust = 2:length(STUDY.cluster)
        STUDY.cluster(clust).topo(find(isnan(STUDY.cluster(clust).topo))) = 0;
    end
    
    % save the clustering results
    all_clusters.(['clusters_' num2str(iteration)]) = STUDY.cluster;
    
end

% save topographies
iterations = 1:max_iterations;

for iteration = iterations
    
    iteration_clusters = all_clusters.(['clusters_' num2str(iteration)]);

    for clust = 3:length(all_clusters.clusters_1)
        
        all_topos(iteration,clust,:) = iteration_clusters(clust).topo(:);
        
    end
end

% correlate

combos = nchoosek(iterations,2);

for combo = 45:size(combos,1)
    start = datetime;
    disp(['Combo ' num2str(combo) '/' num2str(size(combos,1)) ': Cluster solution ' num2str(combos(combo,1)) ' and ' num2str(combos(combo,2))])
    
    for clust = 3:size(all_topos,2)
        % correlate this cluster of this data set with all clusters of the
        % other data set and the other way round at the same time
        for second_clust = 3:size(all_topos,2)
            [there_correlations(clust,second_clust) there_p_values(clust,second_clust)] = corr(squeeze(all_topos(combos(combo,1),clust,:)),squeeze(all_topos(combos(combo,2),second_clust,:)),'type','Pearson');
            [back_correlations(clust,second_clust) back_p_values(clust,second_clust)] = corr(squeeze(all_topos(combos(combo,2),clust,:)),squeeze(all_topos(combos(combo,1),second_clust,:)),'type','Pearson');
        end
        % find the maximum indices showing which topography from one data
        % set correlates best. Goes both ways since the relation is not
        % injective but results may vary. E.g. topographies may be split up in
        % one clustering result. 
        
        % dimensions: 1->index of there max corr, 2->index of back max
        % corr, 3-> there max corr, 4-> back max corr, 5, there max corr p
        % value, 6-> back max corr p value
        index(1) = find(there_correlations(clust,:) == max(there_correlations(clust,:)));
        index(2) = find(back_correlations(clust,:) == max(back_correlations(clust,:)));
        max_correlations(combo,clust,1) = index(1);
        max_correlations(combo,clust,2) = index(2);
        max_correlations(combo,clust,3) = there_correlations(clust,index(1));
        max_correlations(combo,clust,4) = back_correlations(clust,index(2));
        max_correlations(combo,clust,5) = there_p_values(clust,index(1));
        max_correlations(combo,clust,6) = back_p_values(clust,index(2));
    
    end
    thisend = datetime;
    lastduration = thisend - start;
    eta = lastduration * ( size(combos,1) - combo);
    disp('Last duration:')
    disp(lastduration)
    disp('ETA:')
    disp(eta)
end

save('P:\Lukas_Gehrke\studies\Spot_Rotation\data\5_study_level\clusterings_100_dipoles_2_ERSPs_1','clusterings_100_dipoles_2_ERSPs_1','-v7.3')
save('P:\Lukas_Gehrke\studies\Spot_Rotation\data\5_study_level\max_correlations_for_100_clusters','max_correlations','-v7.3')
%% plot stuff

STUDY.cluster = clusterings_100_dipoles_2_ERSPs_1.clusters_1;
std_topoplot(STUDY,ALLEEG,'clusters',3:length(STUDY.cluster));

%% create networks

for combo = 1
    
    firstElement = clusterings_100_dipoles_2_ERSPs_1.(['clusters_' num2str(combos(combo,1))]);
    secondElement = clusterings_100_dipoles_2_ERSPs_1.(['clusters_' num2str(combos(combo,2))]);
    
    firstElement(1:2)=[];
    secondElement(1:2)=[];
    
    this_max_correlations = squeeze(max_correlations(combo,:,[1 2]));
    this_max_correlations(1:2,:) = [];
    
    tupelNetworks = struct();
    tupelNetworksIndex = 1;
    
    while ~isempty(firstElement) || ~isempty(secondElement)
        
        networkOfThisCluster = struct();
        [firstElement, secondElement, this_max_correlations, networkOfThisCluster] = findNetworkForCluster(firstElement,secondElement,this_max_correlations,[1 1],'firstElement',networkOfThisCluster);
        
        tupelNetworks(tupelNetworksIndex).firstElement = networkOfThisCluster.A;
        tupelNetworks(tupelNetworksIndex).secondElement = networkOfThisCluster.B;

        tupelNetworksIndex = tupelNetworksIndex + 1;
        
    end
    
end

function  [newElementA, newElementB, new_max_correlations, newNetworkOfThisCluster] = findNetworkForCluster(firstElement,secondElement,max_correlations,startingCluster,startingElement,networkOfThisCluster)

if strcmp(startingElement,'firstElement')
    
    networkOfThisCluster.A = [networkOfThisCluster.A firstElement(startingCluster)];
    networkOfThisCluster.B = 
    
end

end
