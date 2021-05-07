% pop_bemobil_segment() - Select data bounded by two events. Either keep or
% remove the resulting segment from the data. If multiple segments are
% treated, start_segment and end_segment must be string arrays with the
% same number of elements.
%
% Usage:
%   >>  EEG = bemobil_preprocess( EEG );
%
% Inputs:
%   EEG     - EEGLAB EEG structure
%    
% Outputs:
%   EEG     - EEGLAB EEG structure with boundary events at removed segments
%
% See also: 
%   BEMOBIL_SEGMENT, EEG_FIND_EVENTIDX, POP_SELECT, EEGLAB

function [ EEG ] = pop_bemobil_segment( EEG )

if nargin < 1
	help pop_bemobil_segment;
	return;
end;	

% add preprocess GUI
uigeom_cols = 5;
uigeom_rows = 5;

uigom = {
    {uigeom_cols uigeom_rows [0 0] [1 1]} 
    {uigeom_cols uigeom_rows [1 0] [1 1]}
    {uigeom_cols uigeom_rows [2 0] [1 1]} 
    {uigeom_cols uigeom_rows [3 0] [1 1]}
    {uigeom_cols uigeom_rows [4 0] [1 1]}
    
    {uigeom_cols uigeom_rows [0 1] [1 1]} 
    {uigeom_cols uigeom_rows [1 1] [1 1]}
    {uigeom_cols uigeom_rows [2 1] [1 1]} 
    {uigeom_cols uigeom_rows [3 1] [1 1]}    
    {uigeom_cols uigeom_rows [4 1] [1 1]}    
    
    {uigeom_cols uigeom_rows [0 2] [1 1]} 
    {uigeom_cols uigeom_rows [1 2] [1 1]}
    {uigeom_cols uigeom_rows [2 2] [1 1]} 
    {uigeom_cols uigeom_rows [3 2] [1 1]}
    {uigeom_cols uigeom_rows [4 2] [1 1]}
    
    {uigeom_cols uigeom_rows [0 3] [1 1]} 
    {uigeom_cols uigeom_rows [1 3] [1 1]}
    {uigeom_cols uigeom_rows [2 3] [1 1]} 
    {uigeom_cols uigeom_rows [3 3] [1 1]}
    {uigeom_cols uigeom_rows [4 3] [1 1]}    
    
    {uigeom_cols uigeom_rows [0 4] [1 1]} 
    {uigeom_cols uigeom_rows [1 4] [1 1]}
    {uigeom_cols uigeom_rows [2 4] [1 1]} 
    {uigeom_cols uigeom_rows [3 4] [1 1]}
    {uigeom_cols uigeom_rows [4 4] [1 1]}
    };

uilist = {
    { 'style', 'text', 'string', 'Start Exp Marker:'}
    { 'style', 'edit', 'string', '', 'tag', 'exp_start'}
    { 'style', 'text', 'string', 'End Exp Marker:'}
    { 'style', 'edit', 'string', '', 'tag', 'exp_end'}
    { 'style', 'checkbox', 'string', 'remove', 'value', 0, 'tag', 'keep_or_remove1'}
    
    { 'style', 'text', 'string', 'Start Break 1:'}
    { 'style', 'edit', 'string', '', 'tag', 'break_start1'}
    { 'style', 'text', 'string', 'End Break 1:'}
    { 'style', 'edit', 'string', '', 'tag', 'break_end1'}
    { 'style', 'checkbox', 'string', 'remove', 'value', 1, 'tag', 'keep_or_remove2'}
    
    { 'style', 'text', 'string', 'Start Break 2:'}
    { 'style', 'edit', 'string', '', 'tag', 'break_start2'}
    { 'style', 'text', 'string', 'End Break 2:'}
    { 'style', 'edit', 'string', '', 'tag', 'break_end2'}
    { 'style', 'checkbox', 'string', 'remove', 'value', 1, 'tag', 'keep_or_remove3'}
    
    { 'style', 'text', 'string', 'Start Break 3:'}
    { 'style', 'edit', 'string', '', 'tag', 'break_start3'}
    { 'style', 'text', 'string', 'End Break 3:'}
    { 'style', 'edit', 'string', '', 'tag', 'break_end3'}
    { 'style', 'checkbox', 'string', 'remove', 'value', 1, 'tag', 'keep_or_remove4'}
    
    { 'style', 'text', 'string', 'Start Break 4:'}
    { 'style', 'edit', 'string', '', 'tag', 'break_start4'}
    { 'style', 'text', 'string', 'End Break 4:'}
    { 'style', 'edit', 'string', '', 'tag', 'break_end4'}
    { 'style', 'checkbox', 'string', 'remove', 'value', 1, 'tag', 'keep_or_remove5'}
    };

[ tmp1 tmp2 strhalt structout ] = inputgui('geom', uigom, 'uilist', uilist, ...
    'title', 'Remove or keep segments of data', 'helpcom', 'pophelp(''pop_bemobil_remove_segments'');');
    
if isempty(strhalt)
    disp('Exiting...')
    return; 
end;

% close menu figure after OK
close(gcbf);

break_start = {structout.break_start1 structout.break_start2 structout.break_start3 structout.break_start4};
break_start = break_start(~cellfun('isempty', break_start));

break_end = {structout.break_end1 structout.break_end2 structout.break_end3 structout.break_end4};
break_end = break_end(~cellfun('isempty', break_end));

%if structout.keep_or_remove1 if want to build modular function later

if ~isempty(structout.exp_start)
    EEG = bemobil_segment(EEG, 'keep', structout.exp_start, structout.exp_end);
end

if ~isempty(break_start)
    for i = 1:length(break_start)
        EEG = bemobil_segment(EEG, 'remove', break_start{i}, break_end{i});
    end
end
end