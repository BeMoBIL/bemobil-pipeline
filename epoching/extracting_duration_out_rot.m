%% Extract outward duration

expression_start = 'start_outward_rotation';
expression_end = 'end_outward_rotation'; %all test trials

start_marker = cell(1);
start_marker_time = cell(1);

end_marker = cell(1);
end_marker_time = cell(1);

cnt1 = 1;
cnt2 = 1;

for i = 1:length(EEG.event)
    
    if strcmp(EEG.event(i).type, expression_start)
        start_marker{cnt1} = EEG.event(i).type;
        start_marker_time{cnt1} = EEG.event(i).latency;
        
        cnt1 = cnt1 + 1;
    end    
    
    if strcmp(EEG.event(i).type, expression_end)
        end_marker{cnt2} = EEG.event(i).type;
        end_marker_time{cnt2} = EEG.event(i).latency;
        
        cnt2 = cnt2 + 1;
    end
end

duration = cellfun(@minus, end_marker_time, start_marker_time);
disp(max(duration));

