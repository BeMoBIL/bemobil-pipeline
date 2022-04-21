% bemobil_avref() - Computes the average reference. If EOG channels are declared in channel types, these are excluded.
% If reference channel exists and is declared in EEG.chanlocs(1).ref regular average reference is computed. If no
% reference channel is declared, the full rank average reference is used (see fullrankaveref by Makoto Miakoshi (2017)).
% In both cases the data rank is the number of original channels, so no subtraction of 1 rank for AMICA is necessary!
%
% Usage:
%   >>  [ALLEEG, EEG, CURRENTSET] = bemobil_avref( EEG , ALLEEG, CURRENTSET)
%   >>  [ALLEEG, EEG, CURRENTSET] = bemobil_avref( EEG , ALLEEG, CURRENTSET, out_filename, out_filepath)
%
% Inputs:
%   ALLEEG                  - complete EEGLAB data set structure
%   EEG                     - current EEGLAB EEG structure
%   CURRENTSET              - index of current EEGLAB EEG structure within ALLEEG
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
%   EEGLAB, bemobil_interp, bemobil_avref, pop_interp, pop_reref
% 
% Author: Marius Klug, 2021

function [ALLEEG, EEG, CURRENTSET] = bemobil_avref( EEG , ALLEEG, CURRENTSET, out_filename, out_filepath)

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

%% Compute average reference for all EEG channels

disp('Computing average reference excluding EOG.')
disp('This does NOT reduce the data rank because either full rank average reference is used or the original reference channel was fed back to the data!')

exclude_channels = find(strcmp({EEG.chanlocs.type},'EOG'));
EEG_channels = 1:EEG.nbchan;
EEG_channels(exclude_channels) = [];

if isempty(EEG.chanlocs(1).ref)
    % no ref was declared during preprocessing, use full rank averef (without needing dependency) Apply average
    % reference after adding initial reference, see fullrankaveref by Makoto Miakoshi (2017). This adds an empty new
    % channel, rereferences, then removes the excess channel again, so the rank is still intact. The reference channel,
    % however, is gone and can't be used for analyses, so in case Cz was used as reference during recording, it is
    % inaccessible.
    EEG.nbchan = EEG.nbchan+1;
    EEG.data(end+1,:) = zeros(1, EEG.pnts);

    if size(EEG.chanlocs,1)>1
        EEG.chanlocs = EEG.chanlocs';
    end

    EEG.chanlocs(1,EEG.nbchan).labels = 'initialReference';
    EEG = pop_reref( EEG,[],'keepref','on','exclude',exclude_channels);
    EEG = pop_select( EEG,'nochannel',{'initialReference'});
    
else
    % ref was declared, keep it as channel. this means we have an extra channel, e.g. 129 instead of 128 electrodes, and
    % the former reference carries information. however, the rank is still the previous rank - 1!
    EEG = pop_reref( EEG,[],'keepref','on','exclude',exclude_channels);
end

EEG.etc.bemobil_reref = EEG_channels;

if ~isempty(exclude_channels)
    EEG.etc.bemobil_reref_exclude = EEG.chanlocs(exclude_channels).labels;
else
    EEG.etc.bemobil_reref_exclude = [];
end

% set old REF channel type to EEG now that it is average referenced
[EEG.chanlocs(find(strcmp({EEG.chanlocs.type},'REF'))).type] = deal('EEG');

disp('Rereferencing done.');
