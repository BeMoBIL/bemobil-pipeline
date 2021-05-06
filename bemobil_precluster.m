% bemobil_precluster() - Preparing the STUDY set for clustering. Stores the info in the STUDY struct.
%
% Usage:
%   >>  [STUDY, ALLEEG] = bemobil_precluster(STUDY, ALLEEG, clustering_weights, freqrange, timewindow)
%   >>  [STUDY, ALLEEG] = bemobil_precluster(STUDY, ALLEEG, clustering_weights, freqrange, timewindow, out_filename, out_filepath)
%
% Inputs:
%   ALLEEG                  - complete EEGLAB data set structure
%   STUDY                   - STUDY data set of EEGLAB, which has to be loaded previously
%   clustering_weights      - STRUCT of weights for preclustering, can contain the following fields
%                               dipoles, scalp, spec, erp, ersp
%                             Fields are not required, and weights of 0 will be ignored. At least one field must be
%                             present.
%   freqrange               - frequency range of the ERSPs and spectra to be taken into account
%   timewindow              - time window of ERSPs and ERPs to be taken into account (in ms), can be empty to contain
%                               the complete epoch
%   out_filename            - output filename (OPTIONAL ARGUMENT)
%   out_filepath            - output filepath (OPTIONAL ARGUMENT - File will only be saved on disk
%       if both a name and a path are provided)
%
% Outputs:
%   ALLEEG                  - complete EEGLAB data set structure
%   EEG                     - current EEGLAB EEG structure
%   STUDY                   - STUDY data set of EEGLAB prepared for clustering
%
%   .study data file of current EEGLAB EEG structure stored on disk (OPTIONALLY)
%
% See also:
%   EEGLAB, std_preclust, bemobil_repeated_clustering_and_evaluation
%
% Authors: Marius Klug, 2021

function [STUDY, ALLEEG, EEG] = bemobil_precluster(STUDY, ALLEEG, EEG, clustering_weights, freqrange, timewindow, out_filename, out_filepath)

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

% perform preclustering
command = '[STUDY ALLEEG] = std_preclust(STUDY, ALLEEG, 1';

if isfield(clustering_weights,'dipoles') && clustering_weights.dipoles ~= 0
    command = [command ',{''dipoles'',''weight'',clustering_weights.dipoles}'];
end
if isfield(clustering_weights,'scalp') && clustering_weights.scalp ~= 0
    command = [command ',{''scalp'',''npca'',10,''weight'',clustering_weights.scalp,''abso'',1}'];
end
if isfield(clustering_weights,'spec') && clustering_weights.spec ~= 0
    command = [command ',{''spec'',''npca'',10,''weight'',clustering_weights.spec,''freqrange'',freqrange }'];
end
if isfield(clustering_weights,'erp') && clustering_weights.erp ~= 0
    command = [command ',{''erp'',''npca'',10,''weight'',clustering_weights.erp,''timewindow'',timewindow,''erpfilter'',''''}'];
end
if isfield(clustering_weights,'ersp') && clustering_weights.ersp ~= 0
    command = [command ',{''ersp'',''npca'',10,''freqrange'',freqrange,''timewindow'',timewindow,''norm'',1,''weight'',clustering_weights.ersp}'];
end

command = [command ');'];

eval(command)

% store essential info in STUDY struct for later reading
STUDY.etc.bemobil.clustering.preclustparams.clustering_weights = clustering_weights;
STUDY.etc.bemobil.clustering.preclustparams.freqrange = freqrange;
STUDY.etc.bemobil.clustering.preclustparams.timewindow = timewindow;

% save on disk
if save_file_on_disk
    [STUDY, EEG] = pop_savestudy( STUDY, EEG, 'filename',out_filename,'filepath',out_filepath);
    disp('...done');
end
