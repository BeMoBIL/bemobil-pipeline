% bemobil_segment() - Select data bounded by two events. Either keep or
% remove the resulting segment from the data. If multiple segments are
% treated, start_segment and end_segment must be string arrays with the
% same number of elements. E.g. if start_segment = ['start1' 'start2'] then
% end_segment must be correspondingly ['end1' 'end2']. The routine then
% extracts segment between 'start1' and 'end1' as well as between 'start2' and 'end2'. 
%
% Inputs:
%   keep_or_remove     - 'keep' or 'remove'
%   start_segment      - 'start_marker' of a segment of interest
%   end_segment        - 'end_marker' of a segment
%    
% Outputs:
%   EEG     - EEGLAB EEG structure
%
% See also: 
%   POP_BEMOBIL_SEGMENT, EEG_FIND_EVENTIDX, POP_SELECT, EEGLAB

function [ EEG ] = bemobil_segment(EEG, keep_or_remove, start_segment, end_segment)

if nargin < 1
	help bemobil_segment;
	return;
end;	

[exp_start_event, stridx, start_sec] = eeg_find_eventidx(EEG, start_segment);
if isempty(start_sec) 
    warning('No start index found, using first data point as start!')
    start_sec = 0;
end

[exp_end_event, stridx, end_sec] = eeg_find_eventidx(EEG, end_segment);
if isempty(end_sec)
    warning('No end index found, using last data point as end!')
    end_sec=EEG.times(end)/1000;
end

% get experiment segment
current_segment_sec = [start_sec end_sec];

% extract and retain whole period of 1 condition
if strcmp(keep_or_remove, 'keep')
    EEG = pop_select(EEG, 'time', current_segment_sec);
elseif strcmp(keep_or_remove, 'remove')
    EEG = pop_select(EEG, 'notime', current_segment_sec);
end
end

