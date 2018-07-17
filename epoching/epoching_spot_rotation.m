%% Define parameters

% factors of interest
experiment_condition = ''; % learn or test
movement = '';
trial_nr = '';
direction = '';
speed = '';
angle = '';

%% Find all occurences of a marker_of_interest

expression = 'start_outward_rotation:.*,test'; %all test trials
cnt = 1; % index of occurences

for i = 1:length(EEG.event)
    
    if regexp(EEG.event(i).type, expression)
        
        current_marker = EEG.event(i).type;
        current_marker_split = strsplit(current_marker, ',');
        
        if strcmp(current_marker_split{1}, 'start_outward_rotation:steamvr with PS CAPTURE')
            movement = 'body';
        else
            movement = 'joy';
        end
        
        experiment_condition = current_marker_split{2};
        trial_nr = current_marker_split{3};
        direction = current_marker_split{4};
        speed = current_marker_split{5};
        angle = current_marker_split{6};       
        
        EEG.event(i).experiment_condition = experiment_condition;
        EEG.event(i).movement = movement;
        EEG.event(i).trial_nr = str2double(trial_nr);
        EEG.event(i).direction = direction;
        EEG.event(i).speed = speed;
        EEG.event(i).angle = str2double(angle);
        
    end

    % similar to above extract other markers of interest
    % if regexp()
end

%% Epoch to a marker of interest

% EEG = pop_epoch(EEG);
