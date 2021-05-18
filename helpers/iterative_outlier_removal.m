% ITERATIVE_OUTLIER_REMOVAL - Removal of outliers in a vector based on an iterative sigma threshold approach. Only finds
% the number of removals and the threshold.
% 
% Usage:
%   >>  [n_remove, threshold] = iterative_outlier_removal(data_vector,sd_level,sd_level_increase);
%   
% Required Inputs:
%   data_vector                 - input vector 
% 
% Optional Parameters:
%   sd_level                    - initial sd level for removal (default = 3)
%   sd_level_increase           - sd level increase per iteration (default = 0.1)
% 
% Outputs:
%   n_remove                    - number of removed datapoints
%   threshold                   - final threshold that was used 
%
% See also:
%   nt_zapline_plus, clean_data_with_zapline
%
% Author: Marius Klug, 2021
function [n_remove, threshold] = iterative_outlier_removal(data_vector,sd_level,sd_level_increase)

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
