% bemobil_precluster() - Preparing the STUDY set for clustering. Performs a dimensionality reduction
% to 10 final dimensions using PCA.
%
% Usage:
%   >>  [STUDY, ALLEEG] = bemobil_precluster(STUDY, ALLEEG, clustering_weights, freqrange, timewindow)
%   >>  [STUDY, ALLEEG] = bemobil_precluster(STUDY, ALLEEG, clustering_weights, freqrange, timewindow, out_filename, out_filepath)
%
% Inputs:
%   ALLEEG                  - complete EEGLAB data set structure
%   STUDY                   - STUDY data set of EEGLAB, which has to be loaded previously
%   clustering_weights      - STRUCT of weights for preclustering, containing the fields:
%                               dipoles, scalp_topographies, spectra, ERSPs
%   freqrange               - frequency range of the ERSPs and spectra to be taken into account
%   timewindow              - time window of ERSPs to be taken into account (in ms)
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
%   EEGLAB, std_preclust
%
% Authors: Marius Klug, 2018

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
[STUDY, ALLEEG] = std_preclust(STUDY, ALLEEG, 1,{'scalp' 'npca' 10 'norm' 1 'weight'...
    clustering_weights.scalp_topographies 'abso' 1},{'spec' 'npca' 10 'norm' 1 'weight'...
    clustering_weights.spectra 'freqrange' freqrange },{'dipoles' 'norm' 1 'weight'...
    clustering_weights.dipoles},{'ersp' 'npca' 10 'freqrange' freqrange 'timewindow'...
    timewindow  'norm' 1 'weight' clustering_weights.ERSPs},{'finaldim' 'npca' 10});

% store essential info in STUDY struct for later reading
STUDY.bemobil.clustering.preclustparams.clustering_weights = clustering_weights;
STUDY.bemobil.clustering.preclustparams.freqrange = freqrange;
STUDY.bemobil.clustering.preclustparams.timewindow = timewindow;

% save on disk
if save_file_on_disk
    [STUDY, EEG] = pop_savestudy( STUDY, EEG, 'filename',out_filename,'filepath',out_filepath);
    disp('...done');
end
