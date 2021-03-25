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

% To Do  :  json file names 
%           test files with eloc 
%           test multi run files 
%           make it possible to import only selected data 
%           implement natsortorder

% input check and default value assignment 
%--------------------------------------------------------------------------
if ~isfield(bemobil_config, 'bids_data_folder')
    bemobil_config.bids_data_folder = '1_BIDS-data\';
    warning(['Config field "bids_data_folder" has not been specified- using default folder name ' bemobil_config.bids_data_folder])
end

if ~isfield(bemobil_config, 'other_data_types')
    bemobil_config.other_data_types = {'motion'};
    warning(['Config field "other_data_types" has not been specified- using default value ' bemobil_config.other_data_types{1}])
end

if ~isfield(bemobil_config, 'merge_all_sessions')
    bemobil_config.merge_all_sessions = true;
    warning('Config field "merge_all_sessions" has not been specified- using default value "true"')
end

if ~isfield(bemobil_config, 'resample_freq')
    bemobil_config.resample_freq = 250;
    warning('Config field "resample_freq" has not been specified- using default value 250')
end

% write down the names of other data types then eeg (use BIDS suffix)
otherDataTypes  = bemobil_config.other_data_types; 

% should all sessions be merged? 
mergeAllSessions = bemobil_config.merge_all_sessions; 

% construct the bids data directory 
bidsDir         = fullfile(bemobil_config.study_folder, bemobil_config.bids_data_folder);

% all runs and sessions are merged by default - can be optional e.g., bemobil_config.bids_mergeruns = 1; bemobil_config.bids_mergeses  = 1;
targetDir       = fullfile(bemobil_config.study_folder, bemobil_config.raw_EEGLAB_data_folder);                    % construct using existing config fields
tempDir         = fullfile(targetDir, 'temp_bids'); 

% Import data set saved in BIDS, using a modified version of eeglab plugin 
%--------------------------------------------------------------------------
pop_importbids(bidsDir,'datatypes',otherDataTypes,'outputdir', tempDir);

% Restructure and rename the output of the import function
%--------------------------------------------------------------------------
% list all files and folders in the target folder
subDirList      = dir(tempDir);

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
    
    sesDirList      = dir([tempDir subjectDir]);
    
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
            sesEEGDir       = fullfile(tempDir, subjectDir, sesDir, 'eeg'); 
            sesBEHDir       = fullfile(tempDir, subjectDir, sesDir, 'beh'); 
            sesFilesEEG     = dir(sesEEGDir);
            sesFilesBEH     = dir(sesBEHDir);
            allFiles        = [allFiles sesFilesEEG sesFilesBEH];
        end
        
    else
        
        % for unisession, simply find all files in EEG folder
        eegDir          = fullfile(tempDir, subjectDir, 'eeg'); 
        behDir          = fullfile(tempDir, subjectDir, 'beh'); 
        allFiles        = [dir(eegDir) dir(behDir)] ;
        
    end
    
    % select only .set files
    allFiles = allFiles(contains({allFiles(:).name},'.set')) ;
    
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
        
        if find(strncmp(bidsNameSplit, 'run',3))
            isMultiRun          = true;
            runIndex            = bidsNameSplit{strncmp(bidsNameSplit, 'run',3)}(5:end);
        else
            isMultiRun          = false; 
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
                
                dataDir         = fullfile(tempDir, subjectDir, ['ses-' sessionName] , bidsFolderName);
                
                % move files and then remove the empty eeg folder
                newDir          = fullfile(targetDir, [bemobil_config.filename_prefix num2str(subjectNr)]); 
                
                if ~isdir(newDir)
                    mkdir(newDir)
                end
                
                if isMultiRun
                    bemobilName     = [bemobil_config.filename_prefix, num2str(subjectNr), '_' sessionName '_' bemobilModality, '_rec', runIndex, extension];
                else
                    bemobilName     = [bemobil_config.filename_prefix, num2str(subjectNr), '_' sessionName '_' bemobilModality, extension];
                end
                
                % move the file
                movefile(fullfile(dataDir, bidsName), fullfile(newDir, bemobilName));
     
        else
            dataDir         = fullfile(tempDir, subjectDir, bidsFolderName);
                
            % move files and then remove the empty eeg folder
            newDir          = fullfile(targetDir, [bemobil_config.filename_prefix num2str(subjectNr)]);
            if ~isfolder(newDir)
                mkdir(newDir)
            end
            
            % construct the new file name 
            if isMultiRun
                bemobilName     = [bemobil_config.filename_prefix, num2str(subjectNr), '_' bemobil_config.filenames{1} '_' bemobilModality, '_rec', runIndex '_old', extension];
            else
                bemobilName     = [bemobil_config.filename_prefix, num2str(subjectNr), '_' bemobil_config.filenames{1} '_' bemobilModality '_old', extension];
            end
           
            data    = pop_loadset('filepath', dataDir, 'filename', bidsName); 
            pop_saveset(data, 'filepath', newDir, 'filename', bemobilName); 
            
        end
    end 
end

% delete the temporary directory
disp(['removing ' tempDir])
rmdir(tempDir, 's')

% now synchronize and merge the streams
%--------------------------------------------------------------------------
% list all subject directories
subDirList      = dir(targetDir); 
nameArray       = {subDirList.name};
nameFlagArray   = ~contains(nameArray, '.'); % this is to exclude . and .. folders
subDirList      = subDirList(nameFlagArray);

% first merge all the files if needed
for iSub = 1:numel(subDirList)
    
    % list all files in the subject folder
    subjectFiles = dir(fullfile(targetDir, subDirList(iSub).name));
    
    % iterate over sessions
    for iSes = 1:numel(bemobil_config.filenames)
       
        % find all EEG data
        eegFiles = {subjectFiles(contains({subjectFiles.name}, [bemobil_config.filenames{iSes} '_EEG']) & contains({subjectFiles.name}, '.set')).name};
        eegFiles = sort(eegFiles); % not using natsortorder here - potentially problematic for more than 10 runs? (implausible)  
        
        if numel(eegFiles) > 1
            % if there are multiple runs, merge them all
            EEGMerged = pop_loadset(fullfile(targetDir, subDirList(iSub).name, eegFiles{1}));
            for iFile = 2:numel(eegFiles)
                [EEG2]      = pop_loadset(fullfile(targetDir, subDirList(iSub).name, eegFiles{iFile}));
                EEGMerged   = pop_mergeset(EEGMerged, EEG2);
            end
            EEG                 = EEGMerged;
            EEGFileNameWithRun  = eegFiles{iFile}; 
            nameSplit           = regexp(EEGFileNameWithRun,'_', 'split'); % remove _rec entity
            nameJoined          = join(nameSplit(1:end-1),'_');
            EEGSessionFileName  = [nameJoined{1} '.set'];
        elseif numel(eegFiles) == 1
            EEG                 = pop_loadset(fullfile(targetDir, subDirList(iSub).name, eegFiles{1}));
            EEGSessionFileName  = eegFiles{1};
        else
            warning(['No EEG file found in subject dir ' subDirList(iSub).name ', session ' bemobil_config.filenames{iSes}] )
        end
        
        % resample EEG only when the original sampling rate deviates much
        % from the new
        if abs(bemobil_config.resample_freq - EEG.srate) > 0.01 % this value is not enough
            EEG = pop_resample(EEG, bemobil_config.resample_freq);
        end
        
        EEG.setname = EEG.filename(1:end-8); 
        
        % checkset 
        EEG = eeg_checkset(EEG); 
        
        % save merged EEG file for the session
        EEG = pop_saveset(EEG, 'filename',[EEGSessionFileName(1:end-8) EEGSessionFileName(end-3:end)],'filepath',fullfile(targetDir, subDirList(iSub).name));
        disp(['Saved session file ' EEGSessionFileName(1:end-8) EEGSessionFileName(end-3:end)])
        
        % remove the old EEG gile 
        delete(fullfile(targetDir, subDirList(iSub).name, EEGSessionFileName))
        
        % now iterate over other data types
        for iType = 1:numel(otherDataTypes)
            
            switch otherDataTypes{iType}
                case 'motion'
                    bemobilModality = 'MOCAP';
                otherwise
                    bemobilModality = otherDataTypes{iType};
                    disp(['Unknown modality' otherDataTypes{iType} ', looking for ' otherDataTypes{iType} '.set'])
            end
            
            % find all data of the type
            dataFiles = {subjectFiles(contains({subjectFiles.name}, [bemobil_config.filenames{iSes} '_' bemobilModality]) & contains({subjectFiles.name}, '.set')).name};
            
            if numel(dataFiles) > 1
                % if there are multiple runs, merge them all 
                DATAMerged =  pop_loadset(fullfile(targetDir, subDirList(iSub).name, dataFiles{1}));
                for iFile = 2:numel(dataFiles)
                    DATA2           = pop_loadset(fullfile(targetDir, subDirList(iSub).name, dataFiles{iFile}));
                    DATAMerged      = pop_mergeset(DATAMerged, DATA2);
                end
                DATA             = DATAMerged;
                DATAFileNameWithRun  = dataFiles{iFile};
                nameSplit           = regexp(DATAFileNameWithRun,'_', 'split'); % remove _rec entity
                nameJoined          = join(nameSplit(1:end-1),'_');
                DATASessionFileName  = [nameJoined{1} '.set'];
            elseif numel(dataFiles) == 1
                DATA             = pop_loadset(fullfile(targetDir, subDirList(iSub).name, dataFiles{1}));
                DATASessionFileName  = dataFiles{1};
            else
                warning(['No file of modality ' bemobilModality ' found in subject dir ' subDirList(iSub).name ', session ' bemobil_config.filenames{iSes}] )
            end 
        end
        
        % resample DATA
        if abs(bemobil_config.resample_freq - DATA.srate) > 0.01
            DATA    = pop_resample(DATA, bemobil_config.resample_freq);
        end
        
        % copy events from EEG
        DATA.event = EEG.event;
       
        % checkset 
        DATA = eeg_checkset(DATA); 
        
        % save merged EEG file for the session
        DATA = pop_saveset(DATA, 'filename',[DATASessionFileName(1:end-8) DATASessionFileName(end-3:end)],'filepath',fullfile(targetDir, subDirList(iSub).name));
        disp(['Saved session file ' DATASessionFileName(1:end-8) DATASessionFileName(end-3:end)])
        
        % remove the old .set file 
        delete(fullfile(targetDir, subDirList(iSub).name, DATASessionFileName))
    end
    
    if mergeAllSessions
        
        % merge all EEG data over sessions 
        
        % merge all other data over sessions
    end 
    
end


end
