% bemobil_merge() - Merge multiple .set files and save the created file in
% pre-defined folder, if unspecified, in folder defined by EEG.filepath
%
% Usage:
%   >>  [ ALLEEG EEG CURRENTSET ] = bemobil_merge( ALLEEG, EEG, CURRENTSET );
%   >>  [ ALLEEG EEG CURRENTSET ] = bemobil_merge( ALLEEG, EEG, CURRENTSET, indices );
%   >>  [ ALLEEG EEG CURRENTSET ] = bemobil_merge( ALLEEG, EEG, CURRENTSET, indices, out_filename, out_filepath );
%
% Inputs:
%   ALLEEG                  - complete EEGLAB data set structure
%   EEG                     - current EEGLAB EEG structure
%   CURRENTSET              - index of current EEGLAB EEG structure within ALLEEG
%   indices                 - indices of the sets to merge (OPTIONAL ARGUMENT - GUI WILL OPEN IF NOT PROVIDED)
%   out_filename            - output filename (OPTIONAL ARGUMENT)
%   out_filepath            - output filepath (OPTIONAL ARGUMENT - File will only be saved on disk
%       if both a name and a path are provided)
%
% Outputs:
%   ALLEEG                  - complete EEGLAB data set structure
%   EEG                     - current EEGLAB EEG structure
%   Currentset              - index of current EEGLAB EEG structure within ALLEEG
%
%   .set data file of current EEGLAB EEG structure stored on disk (OPTIONALLY)
%
% See also:
%   EEGLAB, pop_mergeset
% 
% Authors: Lukas Gehrke, Marius Klug, 2017

function [ ALLEEG EEG CURRENTSET ] = bemobil_merge( ALLEEG, EEG, CURRENTSET, indices, out_filename, out_filepath)

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

if exist('indices', 'var')
    EEG = pop_mergeset( ALLEEG, indices );
else
    EEG = pop_mergeset( ALLEEG );
end


[ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, CURRENTSET, 'gui', 'off');

% if data files have been specified by indices, note the file names in EEG.etc
if exist('indices', 'var')
    appended_files={ALLEEG(indices(1)).filename};
    for index=indices(2:end)
        appended_files(index)={ALLEEG(index).filename};
    end
    EEG.etc.appended_files = appended_files;
end

% new data set in EEGLAB
[ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, CURRENTSET, 'gui', 'off');
EEG = eeg_checkset( EEG );

% save on disk
if save_file_on_disk
    EEG = pop_saveset( EEG, 'filename',out_filename,'filepath', out_filepath);
    disp('...done');
end

[ALLEEG EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);