function [corrected_rank_reduction_of_bridges,EEG] = bemobil_find_gel_bridges(EEG,threshold)

if ~exist('threshold','var')
    threshold = 0.99;
end

distances_of_bridges = [];

chancorrs = corr(EEG.data');
triang = triu(chancorrs);
% triang(triang == diag(triang)) = 0;
triang = triang.*~eye(size(triang));
bridges = triang>threshold;
rank_reduction_of_bridges = sum(bridges(:));
[channels_1,channels_2] = ind2sub(size(triang),find(triang>threshold));
bridge_pairs=[channels_1 channels_2];


for i_pair = 1:size(bridge_pairs,1)
    locations_1 = [EEG.chanlocs(bridge_pairs(i_pair,1)).X EEG.chanlocs(bridge_pairs(i_pair,1)).Y EEG.chanlocs(bridge_pairs(i_pair,1)).Z];
    locations_2 = [EEG.chanlocs(bridge_pairs(i_pair,2)).X EEG.chanlocs(bridge_pairs(i_pair,2)).Y EEG.chanlocs(bridge_pairs(i_pair,2)).Z];
    
    distances_of_bridges(i_pair) = sqrt(sum((locations_1 - locations_2).^2));
end


warning(['Found ' num2str(rank_reduction_of_bridges) ' electrode bridges!'])
distances_of_bridges

disp('Only taking into account bridges with a distance less than 35mm!')
corrected_rank_reduction_of_bridges = sum(distances_of_bridges < 35)

EEG.etc.rank_reduction_of_bridges = rank_reduction_of_bridges;
EEG.etc.distances_of_bridges = distances_of_bridges;
EEG.etc.bridge_pairs = bridge_pairs;
EEG.etc.corrected_rank_reduction_of_bridges = corrected_rank_reduction_of_bridges;