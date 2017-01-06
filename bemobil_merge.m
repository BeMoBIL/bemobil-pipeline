% bemobil_merge() - Merge multiple .set files and save the ceated file in
% the same folder as 'merged.set'
%
% Usage:
%   >>  [ ALLEEG EEG CURRENTSET ] = bemobil_merge( ALLEEG, EEG, CURRENTSET );
%   >>  [ ALLEEG EEG CURRENTSET ] = bemobil_merge( ALLEEG, EEG, CURRENTSET, indices );
%
% Inputs:
%   indices      - indices of the sets to merge (OPTIONAL ARGUMENT - GUI WILL OPEN IF NOT PROVIDED)
%
% Outputs:
%   .set file    - EEGLAB compatible .set file
%
% See also:
%   EEGLAB

function [ ALLEEG EEG CURRENTSET ] = bemobil_merge( ALLEEG, EEG, CURRENTSET, indices )

% check if merged file already exist and break if it does
out_filename = 'merged.set';
dir_files = dir(EEG.filepath);
if ismember(out_filename, {dir_files.name})
    error(['Warning: merged file already exists in: ' EEG.filepath '. ' 'Exiting...']);
    %return; use only if warning is provided only on console with disp
end

if exist('indices', 'var')
    EEG = pop_mergeset( ALLEEG, indices );
else
    EEG = pop_mergeset( ALLEEG );
end

% saving data set
[ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, CURRENTSET, 'gui', 'off');
EEG = eeg_checkset( EEG );
EEG = pop_saveset( EEG, 'filename',out_filename,'filepath', [ ALLEEG(CURRENTSET-1).filepath '\']);
disp('...done');
[ALLEEG EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);


end