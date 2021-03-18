function bemobil_bids2set(bemobil_config)
% This function reads in BIDS datasets using the eeglab plugin 
% "bids-matlab-tools" and reorganizes the output to be compatible with 
% BeMoBIL pipeline. For now only EEG data are read and restructured
% To be added :
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

% input check and default value assignment 
%--------------------------------------------------------------------------
if ~isfield(bemobil_config, 'bids_data_folder')
    bemobil_config.bids_data_folder = '1_BIDS-data\';
    warning(['Config field "bids_data_folder" has not been specified- using default folder name ' bemobil_config.bids_data_folder])
end

bidsDir         = fullfile(bemobil_config.study_folder, bemobil_config.bids_data_folder);

% all runs and sessions are merged by default - can be optional e.g., bemobil_config.bids_mergeruns = 1; bemobil_config.bids_mergeses  = 1;
targetDir       = fullfile(bemobil_config.study_folder, bemobil_config.raw_EEGLAB_data_folder);                    % construct using existing config fields

% Import data set saved in BIDS, using the standard eeglab plugin (only EEG)
%--------------------------------------------------------------------------
pop_importbids(bidsDir,'datatypes',{'motion'},'outputdir', targetDir);


% Restructure and rename the output of the import function
%--------------------------------------------------------------------------

% list all files and folders in the target folder
subDirList      = dir(targetDir);

% find all subject folders
dirFlagArray    = [subDirList.isdir];
nameArray       = {subDirList.name};
nameFlagArray   = ~contains(nameArray, '.'); % this is to exclude . and .. folders
subDirList      = subDirList(dirFlagArray & nameFlagArray);

% ToDO : add .json in event.json filename and then remove it - this is to
% let it pass through the bids import function 
% Also, study creation does not work with the eeglab version we are using      

% iterate over all subjects
for iSub = 1:numel(subDirList)
    
    subjectDir      = subDirList(iSub).name;
    
    sesDirList      = dir([targetDir subjectDir]);
    
    % check if data set contains multiple sessions
    isMultiSession = any(contains({sesDirList(:).name},'ses-'));
    
    if isMultiSession
        
        % if multisession, iterate over sessions and concatenate files in EEG folder
        
        dirFlagArray    = [sesDirList.isdir];
        nameArray       = {sesDirList.name};
        nameFlagArray   = ~contains(nameArray, '.'); % this is to exclude . and .. folders
        sesDirList      = sesDirList(dirFlagArray & nameFlagArray);
        
        allFiles        = [];
        for iSes = 1:numel(sesDirList)
            sesDir          = sesDirList(iSes);
            sesEEGDir       = fullfile(targetDir, subjectDir, sesDir, 'eeg'); 
            sesBEHDir       = fullfile(targetDir, subjectDir, sesDir, 'beh'); 
            sesFilesEEG     = dir(sesEEGDir);
            sesFilesBEH     = dir(sesBEHDir);
            allFiles        = [allFiles sesFilesEEG sesFilesBEH];
        end
        
    else
        
        % for unisession, simply find all files in EEG folder
        eegDir          = fullfile(targetDir, subjectDir, 'eeg'); 
        behDir          = fullfile(targetDir, subjectDir, 'beh'); 
        allFiles        = [dir(eegDir) dir(behDir)] ;
        
    end
    
    % select only .set and .fdt files
    allFiles = allFiles(contains({allFiles(:).name},'.set')| contains({allFiles(:).name},'.fdt')) ;
    
    for iFile = 1:numel(allFiles)
        
        % rename files to bemobil convention (only eeg files for now)
        bidsName        = allFiles(iFile).name;                             % 'sub-003_task-VirtualNavigation_eeg.set';
        bidsNameSplit   = regexp(bidsName, '_', 'split');
        subjectNr       = str2double(bidsNameSplit{1}(5:end));
        bidsModality    = bidsNameSplit{end}(1:end-4);                        % this string includes modality and extension
        extension       = bidsNameSplit{end}(end-3:end);
        
        if find(strncmp(bidsNameSplit, 'ses',3))
            isMultiSession      = 1; 
            sessionName         = bidsNameSplit{strncmp(bidsNameSplit, 'ses',3)}(5:end);
        end
        
        switch bidsModality
            case 'eeg'
                bemobilModality = upper(bidsModality);                      % use string 'EEG' for eeg data
                bidsFolderName  = 'eeg'; 
            case 'motion'
                bemobilModality = 'MOCAP';
                bidsFolderName = 'beh'; 
            otherwise
                bemobilModality = bidsModality;
                disp(['Unknown modality' bidsModality ' saved as ' bidsModality '.set'])
        end
        
        if isMultiSession
                
                dataDir         = fullfile(targetDir, subjectDir, ['ses-' sessionName] , bidsFolderName);
                
                % move files and then remove the empty eeg folder
                newDir     = fullfile(targetDir, [bemobil_config.filename_prefix num2str(subjectNr)]); 
                if ~isdir(newDir)
                    mkdir(newDir)
                end
                
                bemobilName     = [bemobil_config.filename_prefix num2str(subjectNr), '_' sessionName '_' bemobilModality extension];
    
                % move the file
                movefile(fullfile(dataDir, bidsName), fullfile(newDir, bemobilName));
                if numel(dir(dataDir)) == 2
                    rmdir(dataDir)
                end
             
            sesDir =  fullfile(targetDir, subjectDir, ['ses-' sessionName]);
            if numel(dir(sesDir)) == 2
                rmdir(sesDir)
            end
        else
            dataDir         = fullfile(targetDir, subjectDir, bidsFolderName);
                
            % move files and then remove the empty eeg folder
            newDir     = fullfile(targetDir, [bemobil_config.filename_prefix num2str(subjectNr)]);
            if ~isfolder(newDir)
                mkdir(newDir)
            end
            % construct the name with filename in the middle 
            bemobilName     = [bemobil_config.filename_prefix num2str(subjectNr), '_' bemobil_config.filenames{1} '_' bemobilModality extension];
            movefile( fullfile(dataDir, bidsName), fullfile(newDir, bemobilName));
            if numel(dir(dataDir)) == 2
                rmdir(dataDir)
            end
        end
    end
    
    % if empty, also remove the subject directory    
    if numel(dir([targetDir, subjectDir])) == 2
        rmdir([targetDir, subjectDir])
    end
end

% now synchronize and merge the streams 
%--------------------------------------------------------------------------
% synchronization has to be done using BIDS srate information 

% find the sampling rate of the eeg files

% find the sampling rate of the matching motion files

% upsample motion data if necesary

% merge files and delete the unnecessary pieces of files  

end
