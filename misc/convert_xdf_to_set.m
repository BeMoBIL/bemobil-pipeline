% convert_xdf_to_set() - Convert an XDF file or a folder containing
% multiple XDF files to EEGLAB .set files
%
% Usage:
%   >> [ALLEEG EEG CURRENTSET] = convert_xdf_to_set(ALLEEG, filename, filepath_input, filepath_output, prefix )
%
% Inputs:
%   filename        - filename of a *.xdf file
%   filepath_input  - filepath to a *.xdf file
%   filepath_output - filepath for the *.set file
%    
% Outputs:
%   .set file    - EEGLAB compatible .set file
%
% See also: 
%   LOAD_XDF, POP_LOAD_XDF, POP_SAVESET, EEGLAB

function [ALLEEG EEG CURRENTSET] = convert_xdf_to_set(ALLEEG, filename, filepath_input, filepath_output, prefix )

if nargin < 1
	help convert_xdf_to_set;
	return;
end;

if ~exist('filepath_output','var')
    filepath_output = filepath_input;
end

if ~exist('prefix','var')
	prefix = [];
end

mkdir(filepath_output)

if isempty(filename)
    % if a single file is not provided, find all *.xdf
    % files in session_path
    data_to_process_tmp = dir([filepath_input '\*.xdf']);
    data_to_process = {data_to_process_tmp.name};
else
    data_to_process = strcat(filename,'.xdf');
end

% check if any data files are ready to be processed
if isempty(data_to_process)
    disp('Warning!!! No data present to be preprocessed. Exiting...');
    return;
end

% 1. load data
for i = 1:length(data_to_process)
    disp(['Now loading dataset: ' num2str(i) ' in session folder']);
    
    setname = strsplit(data_to_process{i},'.');
    setname = strcat(prefix, setname{1});
    EEG = pop_loadxdf(strcat(filepath_input, '/', data_to_process{i}), 'streamtype', 'EEG', 'exclude_markerstreams', {});
    [ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 1,'setname',[setname '_EEG'],'savenew',[filepath_output '\' setname '_EEG.set'],'gui','off'); 

end
end

