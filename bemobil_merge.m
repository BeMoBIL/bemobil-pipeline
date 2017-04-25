% bemobil_merge() - Merge multiple .set files and save the ceated file in
% the same folder as 'merged.set'
%
% Usage:
%   >>  [ ALLEEG EEG CURRENTSET ] = bemobil_merge( ALLEEG, EEG, CURRENTSET );
%   >>  [ ALLEEG EEG CURRENTSET ] = bemobil_merge( ALLEEG, EEG, CURRENTSET, indices );
%
% Inputs:
%   indices      - indices of the sets to merge (OPTIONAL ARGUMENT - GUI WILL OPEN IF NOT PROVIDED)
%   out_filename - output filename for the merged file (OPTIONAL ARGUMENT - DEFAULT VALUE IS 'merged.set')
%   out_filepath - output filepath for the merged file (OPTIONAL ARGUMENT - DEFAULT VALUE IS SAME FOLDER)
%
% Outputs:
%   EEG          - EEGLAB EEG data of merged datasets
%   .set file    - EEGLAB compatible .set file is saved on disk
%
% See also:
%   EEGLAB
% 
% Authors: Lukas Gehrke, Marius Klug, Friederike Hohlefeld, 2017

function [ ALLEEG EEG CURRENTSET ] = bemobil_merge( ALLEEG, EEG, CURRENTSET, indices, out_filename, out_filepath)

if ~exist('out_filename', 'var') out_filename = 'merged.set'; end;
if ~exist('out_filepath', 'var') out_filepath = EEG.filepath; end;

% check if merged file already exist and break if it does
mkdir(out_filepath); % make sure that folder exists, nothing happens if so
dir_files = dir(out_filepath);
if ismember(out_filename, {dir_files.name})
    error(['Warning: ' out_filename ' file already exists in: ' out_filepath '. ' 'Exiting...']);
    %return; use only if warning is provided only on console with disp
end

if exist('indices', 'var')
    EEG = pop_mergeset( ALLEEG, indices );
else
    EEG = pop_mergeset( ALLEEG );
end


[ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, CURRENTSET, 'gui', 'off');

% if data files have been specified by indices, note the file names in EEG.etc
if exist('indices', 'var')
    appended_files = ALLEEG(indices(1)).filename;
    for index=indices(2:end)
        appended_files = [appended_files ', ' ALLEEG(index).filename];
    end
    EEG.etc.appended_files = appended_files;
end

EEG = eeg_checkset( EEG );
EEG = pop_saveset( EEG, 'filename',out_filename,'filepath', out_filepath);
disp('...done');
[ALLEEG EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);


end