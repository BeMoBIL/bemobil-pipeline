function bemobil_bids2set(bemobil_config)
% This function reads in BIDS datasets using the eeglab plugin 
% "bids-matlab-tools" and reorganizes the output to be compatible with 
% BeMoBIL pipeline. For now only EEG data are read and restructured
% To be added :
%           reading in Motion data or data of other modalities
%           support separate output files for multi-run and multi-session
%
% Usage
%       bemobil_bids2set(bemobil_config)
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
%               (link to be provided)
%       bva-io for brain vision data :
%               https://github.com/arnodelorme/bva-io
%
% author : seinjeung@gmail.com
%--------------------------------------------------------------------------

% path to the folder containing BIDS data
bemobil_config.bids_folder                  =  '\\stor1\projects\Sein_Jeung\Project_Virtual_Navigation\Virtual_Navigation_data\data_E1\rawdata';

% this is the target folder in which the data will be written
bemobil_config.study_folder                 = '\\stor1\projects\Sein_Jeung\Project_Virtual_Navigation\Virtual_Navigation_data\data_E1';
bemobil_config.raw_EEGLAB_data_folder       = '2_basic-EEGLAB-test\';

% name to be used for the final merged EEG file
bemobil_config.filenames                    = {'VN_E1'};

% all runs and sessions are merged by default - can be optional e.g., bemobil_config.bids_mergeruns = 1; bemobil_config.bids_mergeses  = 1;
bidsDir         = bemobil_config.bids_folder;
targetDir       = fullfile(bemobil_config.study_folder, bemobil_config.raw_EEGLAB_data_folder);                    % construct using existing config fields

% Import data set in BIDS using the standard eeglab plugin (only EEG)
%--------------------------------------------------------------------------
pop_importbids(bidsDir, 'outputdir', targetDir);

% Restructure and rename the output of the import function
%--------------------------------------------------------------------------

% list all files and folders in the target folder
subDirList      = dir(targetDir);

% find all subject folders
dirFlagArray    = [subDirList.isdir];
nameArray       = {subDirList.name};
nameFlagArray   = ~contains(nameArray, '.'); % this is to exclude . and .. folders
subDirList      = subDirList(dirFlagArray & nameFlagArray);

% iterate over all subjects
for iSub = 1:numel(subDirList)
    
    subjectDir      = subDirList(iSub).name;
    
    % check if data set contains multiple sessions
    isMultiSession = any(strcmp(subjectDir(:).name(1:3),'ses-'));
    
    if isMultiSession
        
        % if multisession, iterate over sessions and concatenate files in EEG folder
        sesDirList      = dir(subjectDir);
        dirFlagArray    = [sesDirList.isdir];
        nameArray       = {sesDirList.name};
        nameFlagArray   = ~contains(nameArray, '.'); % this is to exclude . and .. folders
        sesDirList      = sesDirList(dirFlagArray & nameFlagArray);
        
        eegFiles        = [];
        for iSes = 1:numel(sesDirList)
            sesDir      = sesDirList(iSes);
            sesFiles    = dir(sesDir);
            eegFiles    = [eegFiles sesFiles];
        end
        
    else
        
        % for unisession, simply find all files in EEG folder
        eegFiles       = dir([subjectDir '\eeg']);
        
    end
    
    % select only .set and .fdt files
    eegFiles = eegFiles(strcmp(eegFiles(:).name(end-3:end),'.set')|| strcmp(eegFiles(:).name(end-3:end),'.fdt')) ;
    
    for iFile = 1:numel(eegFiles)
        
        % rename files to bemobil convention (only eeg files for now)
        bidsName        = eegFiles(iFile).name;                             % 'sub-003_task-VirtualNavigation_eeg.set';
        bidsNameSplit   = regexp(bidsName, '_', 'split');
        subjectNr       = str2double(bidsNameSplit{1}(5:end));
        bidsModality    = bidsNameSplit{end}(1:end-4);                      % this string includes modality and extension
        extension       = bidsNameSplit{end}(end-4:end);
        
        switch bidsModality
            case 'eeg'
                bemobilModality = upper(bidsModality);                      % use string 'EEG' for eeg data
            case 'motion'
                disp('Found motion data in .set format - not implemented yet')
            otherwise
                bemobilModality = bidsModality;
                disp(['Unknown modality' bidsModality ' saved as ' bidsModality '.set'])
        end
        
        bemobilName     = [bemobil_config.filename_prefix num2str(subjectNr) '_' strjoin(bidsNameSplit(2:end-1),'_') '_' bemobilModality extension];
        movefile( fullfile(projectdir, bidsName), fullfile(projectdir, bemobilName));
        
    end
end

end
