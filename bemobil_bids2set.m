function bemobil_bids2set(bemobil_config)
% This function reads in data sets that are organized according to BIDS
% and reorganizes the output to be compatible with bemobil pipeline
% To be added : synchronization between streams 
%
% Usage 
%       bemobil_bids2set(config)
% 
% In 
%       config
%       see help bemobil_config documentation 
%
% Out 
%       none
%       reorganizes data on disk 
%
% required plugins
%       modified version of SCCN bids-matlab-tools :
%       (link to be provided)
%       bva-io : https://github.com/arnodelorme/bva-io
%
% author : seinjeung@gmail.com
%--------------------------------------------------------------------------
bemobil_config.study_folder                 = '\\stor1\projects\Sein_Jeung\Project_Virtual_Navigation\Virtual_Navigation_data\data_E1'; 
bemobil_config.filenames                    = {'VN_E1'};
bemobil_config.raw_EEGLAB_data_folder       = '2_basic-EEGLAB-test\';
bemobil_config.bids_folder                  =  '\\stor1\projects\Sein_Jeung\Project_Virtual_Navigation\Virtual_Navigation_data\data_E1\rawdata';
bemobil_config.data_types                   = [{'motion'}]; 

bidsDir         = bemobil_config.bids_folder; % addition to bemobil_config
datatypes       = bemobil_config.data_types; % addition to bemobil_config 
targetDir       = fullfile(bemobil_config.study_folder, bemobil_config.raw_EEGLAB_data_folder);                    % construct using existing config fields

% Import data set in BIDS using the standard eeglab plugin 
%--------------------------------------------------------------------------
pop_importbids(bidsDir, 'datatypes', datatypes, 'outputdir', targetDir); 

% Restructure and rename the output of import function
%--------------------------------------------------------------------------
% iterate over all subjects
for iSub = 1:nSubjects
    
    % 1. rename subject files according to bemobil convention
    subjectDir      = subjectDir(iSub);
    
    % 2. move data files from modality specific folders to subject folders
    bidsFiles       = dir(subjectDir);
    
    % 3. rename files to bemobil convention
    bidsName        = 'sub-003_task-VirtualNavigation_eeg.set';
    bidsNameSplit   = regexp(bidsName, '_', 'split');
    subjectNr       = str2double(bidsNameSplit{1}(5:end));
    bidsModality    = bidsNameSplit{3}(1:end-4);                             % this string includes modality and extension
    extension       = bidsNameSplit{3}(end-4:end);
    taskLabel       = '';
    acqLabel        = '';
    
    switch bidsModality
        case 'eeg'
            bemobilModality = upper(bidsModality);                          % use string 'EEG' for eeg data
        case 'motion'
            bemobilModality = upper(bidsModality);                          % use string 'MOTION' for motion data
        otherwise
            bemobilModality = bidsModality;
            disp(['Unknown modality' bidsModality 'saved as ' bidsModality '.set'])
    end
    
    bemobilName     = [bemobil_config.filename_prefix num2str(subjectNr) '_' taskLabel '_' acqLabel '_' bemobilModality extension];
    movefile( fullfile(projectdir, bidsName), fullfile(projectdir, bemobilName) );
    
end

end
