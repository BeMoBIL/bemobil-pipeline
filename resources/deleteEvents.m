%% deleting events

temp = EEG.event;

% delete all events that are not purely alphabetical
indNumericalChars = cellfun(@(x) sum(isstrprop(x, 'digit'))>1, {temp.type});

temp(indNumericalChars)=[];


% delete all events that are neither start, stop, begin, end,
% movementonset, or boundary
indMovementOnset =  ~cellfun(@isempty,strfind(lower({temp.type}),'movementonset'));
indStarts =  ~cellfun(@isempty,strfind(lower({temp.type}),'start'));
indStops =  ~cellfun(@isempty,strfind(lower({temp.type}),'stop'));
indBegins =  ~cellfun(@isempty,strfind(lower({temp.type}),'begin'));
indEnds =  ~cellfun(@isempty,strfind(lower({temp.type}),'end'));
indBoundaries =  ~cellfun(@isempty,strfind(lower({temp.type}),'boundary'));

indAll = indMovementOnset + indStarts + indStops + indBegins + indEnds + indBoundaries;
temp(~indAll)=[];

EEG.event = temp;
EEG = eeg_checkset( EEG );
eeglab redraw