%eeg_find_eventidx % find index number of event from EEG struct (load with EEGlab)
% string can be cell array or single string

function [events, str_idx, evt_sec] = eeg_find_eventidx(EEG, string)

events = '';
evt_sec = '';

S = {EEG.event.type};
[dummy, str_idx] = ismember(string,S);

if str_idx == 0
    warning(['''' string ''' event not found in dataset'])
    return;
end

events = EEG.event(str_idx);
evt_sec = (cell2mat({events(1:length(events)).latency})-1)/EEG.srate;

end