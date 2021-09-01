% bemobil_create_peak_events() - Inserts events for peak events like eyeblinks or heartbeats, based on the activity of a
% source component of the spatial filter (like Independent Component Analysis). The threshold of a peak is defined as 
% 1/5 of the maximum of the activity in the component. A refractory period after the peak event can be defined, otherwise
% a new peak can be detected as soon as the previous' peak's activity is back below half the threshold. 
%
% Usage:
%   >>  [ ALLEEG EEG CURRENTSET ] = bemobil_create_peak_events(ALLEEG, EEG, CURRENTSET, component_to_use, event_name, refractory_period)
%   >>  [ ALLEEG EEG CURRENTSET ] = bemobil_create_peak_events(ALLEEG, EEG, CURRENTSET, component_to_use, event_name, refractory_period, out_filename, out_filepath)
%
% Inputs:
%   ALLEEG                  - complete EEGLAB data set structure
%   EEG                     - current EEGLAB EEG structure
%   CURRENTSET              - index of current EEGLAB EEG structure within ALLEEG
%   component_to_use        - the source component that should be used for this analysis (integer value)
%   event_name              - the name of the event to be inserted (string value)
%   refractory_period       - minimum time (in ms) between two peaks. Everything in between will not be evaluated. If
%                               empty, the next peak can be detected as soon as the previous' peak's activity is below
%                               half the threshold again
%   out_filename            - output filename (OPTIONAL ARGUMENT)
%   out_filepath            - output filepath (OPTIONAL ARGUMENT - File will only be saved on disk
%       if both a name and a path are provided)
%
% Outputs:
%   ALLEEG                  - complete EEGLAB data set structure including the new events
%   EEG                     - current EEGLAB EEG structure
%   Currentset              - index of current EEGLAB EEG structure within ALLEEG
%
%   .set data file of current EEGLAB EEG structure stored on disk (OPTIONALLY)
%
% See also:
%   EEGLAB, pop_detecteyemovements, detecteyemovements, addevents
% 
% Authors: Marius Klug, 2018

function [ ALLEEG EEG CURRENTSET ] = bemobil_create_peak_events(ALLEEG, EEG, CURRENTSET, component_to_use, event_name, refractory_period, out_filename, out_filepath)

% only save a file on disk if both a name and a path are provided
save_file_on_disk = (exist('out_filename', 'var') && exist('out_filepath', 'var'));

% check if file already exist and show warning if it does
if save_file_on_disk
    mkdir(out_filepath); % make sure that folder exists, nothing happens if so
    dir_files = dir(out_filepath);
    if ismember(out_filename, {dir_files.name})
        warning([out_filename ' file already exists in: ' out_filepath '. File will be overwritten...']);
    end
end

component_activity = EEG.icaact(component_to_use,:);

% sorted = sort(component_activity);
% threshold = sorted(round(length(sorted)*(1-(10/(60*10)))))

% here the threshold is defined -> 1/5 of the maximum in the whole data set
threshold = max(component_activity)*0.20;

event_latencies = [];

timepoint = 1;

% loop through all timepoints
while timepoint<=length(component_activity)
   
    if component_activity(timepoint)>threshold
        % peak detected
        
        starting_timepoint = timepoint;
        
        % find the end of this peak's activity - this is defined as half the threshold to ensure that the activity is
        % actually down to normal and not just fluctuating around the threshold.
        while component_activity(timepoint)>threshold*0.5
            timepoint = timepoint + 1;
        end
        ending_timepoint = timepoint;
        
        maximum_timepoint = starting_timepoint - 1 +...
            find(component_activity(starting_timepoint:ending_timepoint) ==...
            max(component_activity(starting_timepoint:ending_timepoint)));
        
        event_latencies(end+1) = maximum_timepoint;
        
        if ~isempty(refractory_period)
            
            % add minimum diration (in ms) to the maximum timepoint and use that one instead of the old end-of-peak
            % timepoint
            timepoint = maximum_timepoint + round((refractory_period/(1/EEG.srate))/1000);
            
        end
    else
        % no peak detected
        timepoint = timepoint + 1;        
    end
    
end

% make ur-events (necessary for addevents() )
% EEG = eeg_checkset(EEG,'makeur');


for latency = event_latencies
    i = numel(EEG.event) + 1;
    EEG.event(i).type = event_name;
    EEG.event(i).latency = latency;
    EEG.event(i).duration = 1/EEG.srate;
end

EEG = eeg_checkset(EEG, 'eventconsistency');

% % piggyback on a function created in the EYE-EEG toolbox because I didn't want to code this myself
% EEG = addevents(EEG,event_latencies',{'latency'},event_name);

% [ALLEEG EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);pop_eegplot( EEG, 0, 1, 1);
% 
% EEG.event(2:end) = [];

% new data set in EEGLAB
[ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, CURRENTSET, 'gui', 'off');
EEG = eeg_checkset( EEG );

% save on disk
if save_file_on_disk
    EEG = pop_saveset( EEG, 'filename',out_filename,'filepath', out_filepath);
    disp('...done');
end

[ALLEEG EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);
