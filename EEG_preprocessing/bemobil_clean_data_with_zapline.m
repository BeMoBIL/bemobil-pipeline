function [EEG, plothandles] = bemobil_clean_data_with_zapline(EEG, zaplineConfig)
% Small wrapper for zapline-plus to be used with EEGLAB EEG structs, that makes sure all info is stored, and uses a
% struct as input. Removial of frequency artifacts using ZapLine to remove line noise from EEG/MEG data. Adds automatic
% detection of the number of components to remove, and chunks the data into segments to account for nonstationarities.
% Based on: de Cheveigne, A. (2020) ZapLine: a simple and effective method to remove power line artifacts. Neuroimage,
% 1, 1-13.
% 
% Requires noisetools to be installed: http://audition.ens.fr/adc/NoiseTools/
% 
% Usage:
% >> [EEG, plothandles] = bemobil_clean_data_with_zapline(EEG, zaplineConfig)
% 
% Required Inputs:
%   EEG                     - EEGLAB data struct
%   zaplineConfig.linefreqs - vector with one or more line frequencies to be removed
% 
% Optional parameters for the zaplineConfig:
%   adaptiveNremove         - bool. if automatic adaptation of removal should be used. (default = 1)
%   fixedNremove            - numerical vector. fixed number of removed components. if adaptive removal is used, this 
%                               will be the minimum. can be either a scalar (then it is used for all line freqs) or a 
%                               vector of the same length as the linefreqs (then individual fixed n remove will be
%                               used). (default = 0)
%   chunkLength             - numerical. length of chunks to be cleaned in seconds. if set to 0, no chunks will be used. 
%                               (default = 30)
%   plotResults             - bool. if plot should be created. takes time to compute the spectrum. (default = 1)
%   figBase                 - integer. figure number to be created and plotted in. each iteration of linefreqs increases 
%                               this number by 1. (default = 100)
%   nfft                    - numerical. fft window size for computing the spectrum. (default = 512)
%   nkeep                   - integer. PCA reduction of components before removal. (default = round(20+size(data,2)/4))
%   initialSigma            - numerical. initial iterative outlier detection sigma threshold. (default = 3)
%   sigmaIncrease           - numerical. iterative outlier detection sigma threshold increase per iteration (to ensure 
%                               convergence). (default = 0.1)
% 
% Outputs:
%   EEG                     - EEGLAB data struct cleaned with zapline-plus
%   plothandles             - vector of handles to the created figures
%
% 
% Example:
%   [EEG, plothandles] = bemobil_clean_data_with_zapline(EEG, struct('linefreqs',50))
%
% See also:
%   clean_data_with_zapline, nt_zapline_plus, iterative_outlier_removal
%
% Author: Marius Klug, 2021

zaplineFields = fieldnames(zaplineConfig);

zaplineFields(contains(zaplineFields,'linefreqs')) = [];

varargin = {};
for i_fieldname = 1:length(zaplineFields)
    
    varargin{1+(i_fieldname-1)*2} = zaplineFields{i_fieldname};
    varargin{2+(i_fieldname-1)*2} = zaplineConfig.(zaplineFields{i_fieldname});
    
end

% run zapline subfunction with applied parameters
[EEG.data, EEG.etc.zapline.NremoveFinal, EEG.etc.zapline.Scores, EEG.etc.zapline.config, plothandles] =...
    clean_data_with_zapline(EEG.data, EEG.srate, zaplineConfig.linefreqs, varargin{:});
