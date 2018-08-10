if ~isfield(EEG, 'event')
    % creating new empty event structure
    EEG.event = struct('type', {}, 'latency', {}, 'duration', {});
end

% adding peak events
for latency = event_latencies
    i = numel(EEG.event) + 1;
    EEG.event(i).type = 'photo-on';
    EEG.event(i).latency = latency;
    EEG.event(i).duration = 1/EEG.srate;
end

EEG = eeg_checkset(EEG, 'eventconsistency');