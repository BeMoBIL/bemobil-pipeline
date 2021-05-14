% bemobil_split_MoBI_set() - Splits a MoBI dataset into it's unique data types, as stored in EEG.chanlocs.type. The
% mobilab export2eeglab() function stores data types (e.g. EEG, motion, EyeTracking) automatically.
%
% Usage:
%   >>  [ALLEEG, EEG, CURRENTSET, EEG_split_sets] = bemobil_split_MoBI_set(ALLEEG, MoBI_EEG, CURRENTSET);
%
% Inputs:
%   ALLEEG                  - complete EEGLAB data set structure
%   EEG                     - current EEGLAB EEG structure with MoBI data type stored in EEG.chanlocs.type
%   CURRENTSET              - index of current EEGLAB EEG structure within ALLEEG
%
% Outputs:
%   ALLEEG                  - complete EEGLAB data set structure
%   EEG                     - last EEGLAB EEG structure
%   CURRENTSET              - index of last EEGLAB EEG structure within ALLEEG
%   EEG_split_sets          - struct array with the split data sets
%
% See also:
%   EEGLAB, pop_select
%
% Authors: Marius Klug, 2018

function [ALLEEG, EEG, CURRENTSET, EEG_split_sets] = bemobil_split_MoBI_set(ALLEEG, MoBI_EEG, CURRENTSET)

disp('Splitting MoBI dataset in unique types...')

unique_types = unique({MoBI_EEG.chanlocs.type});

for unique_type = 1:length(unique_types)
    
    indices = find(strcmp({MoBI_EEG.chanlocs.type},unique_types(unique_type)));
    
    fprintf('Type "%s": %d of %d channels. ',unique_types{unique_type}, length(indices) , length(MoBI_EEG.chanlocs))
    
    EEG_split_sets(unique_type) = pop_select( MoBI_EEG,'channel',indices);
    
    EEG_split_sets(unique_type) = eeg_checkset( EEG_split_sets(unique_type) );
    
    % new data set in EEGLAB
    [ALLEEG, EEG, CURRENTSET] = pop_newset(ALLEEG, EEG_split_sets(unique_type), CURRENTSET, 'gui', 'off');
    
end

disp('...done!')