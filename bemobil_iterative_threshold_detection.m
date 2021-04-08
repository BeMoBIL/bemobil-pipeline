function [n_remove, threshold] = bemobil_iterative_threshold_detection(data_vector,sd_level,sd_level_increase)

if ~exist('sd_level','var')
	sd_level = 3;
end

if ~exist('sd_level_increase','var')
	sd_level_increase = 0.1;
end

threshold_old = max(data_vector);
threshold = mean(data_vector)+sd_level*std(data_vector);
n_remove = 0;

while threshold < threshold_old
    
	flagged_points = data_vector>threshold;
	
	data_vector(flagged_points) = [];

	n_remove = n_remove + sum(flagged_points);
    
    sd_level = sd_level + sd_level_increase;
    
    threshold_old = threshold;
	threshold = mean(data_vector)+sd_level*std(data_vector);
end
