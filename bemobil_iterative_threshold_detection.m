function [n_remove, threshold] = bemobil_iterative_threshold_detection(data_vector,sd_level)

if ~exist('sd_level','var')
	sd_level = 2.5;
end

sd_level_increase = sd_level*0.075;

threshold_old = max(data_vector);
threshold = mean(data_vector)+sd_level*std(data_vector);
n_remove = 0;

while threshold < threshold_old
	
    threshold_old = threshold;
    
	flagged_points = data_vector>mean(data_vector)+sd_level*std(data_vector);
	
	data_vector(flagged_points) = [];
	
	threshold = mean(data_vector)+3*std(data_vector);

	n_remove = n_remove + sum(flagged_points);
    
    sd_level = sd_level + sd_level_increase;
end