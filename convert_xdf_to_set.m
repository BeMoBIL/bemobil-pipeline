% convert_xdf_to_set() - Convert an XDF file or a folder containing
% multiple XDF files to EEGLAB .set files
%
% Usage:
%   >>  EEG = bemobil_preprocess( filename, filepath );
%
% Inputs:
%   filename     - filename of a *.xdf file
%   filepath     - filepath to a *.xdf file
%    
% Outputs:
%   .set file    - EEGLAB compatible .set file
%
% See also: 
%   LOAD_XDF, POP_LOAD_XDF, POP_SAVESET, EEGLAB

function [ EEG ] = convert_xdf_to_set( filename, filepath )

if nargin < 1
	help convert_xdf_to_set;
	return;
end;

if isempty(filename)
    % if a single file is not provided, find all *.xdf
    % files in session_path
    data_to_process_tmp = dir([filepath '/*.xdf']);
    data_to_process = {data_to_process_tmp.name};
else
    data_to_process = {strcat(filepath, '/', filename)};
end

% check if any data files are ready to be processed
if isempty(data_to_process)
    disp('Warning!!! No data present to be preprocessed. Exiting...');
    return;
end

% 1. load data
for i = 1:length(data_to_process)
    disp(['Now loading dataset: ' num2str(i) ' in session folder']);

    if length(data_to_process) > 1
        EEG = pop_loadxdf(strcat(filepath, '/', data_to_process{i}), 'streamtype', 'EEG', 'exclude_markerstreams', {});
        pop_saveset(EEG, data_to_process{i}, filepath);
    else
        EEG = pop_loadxdf(data_to_process{i}, 'streamtype', 'EEG', 'exclude_markerstreams', {});
        pop_saveset(EEG, filename, filepath);
    end
end
end

