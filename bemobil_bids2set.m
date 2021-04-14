function bemobil_bids2set(bemobil_config, numericalIDs)
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

% add natsortfiles to path
[filepath,~,~] = fileparts(which('bemobil_bids2set')); 
addpath(fullfile(filepath, 'resources', 'natsortfiles'))


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

% check if final session file names are there - if so, skip
% subjectNr = str2double(subDirList(iSub).name(numel(bemobil_config.filename_prefix) + 1:end));
% sessionNameEEG = [bemobil_config.filename_prefix, num2str(subjectNr), '_', bemobil_config.filenames{iSes} '_EEG.set'];
% 
% if isfile(fullfile(targetDir, subDirList(iSub).name, sessionNameEEG))
%     disp(['File ' sessionNameEEG ' already found - skipping import for this participant'])
% end

% Import data set saved in BIDS, using a modified version of eeglab plugin 
%--------------------------------------------------------------------------
pop_importbids_mobi(bidsDir,'datatypes',otherDataTypes,'outputdir', tempDir, 'participants', numericalIDs);

% Restructure and rename the output of the import function
%--------------------------------------------------------------------------
% list all files and folders in the target folder
subDirList      = dir(tempDir);

% find all subject folders
dirFlagArray    = [subDirList.isdir];
nameArray       = {subDirList.name};
nameFlagArray   = ~contains(nameArray, '.'); % this is to exclude . and .. folders
subDirList      = subDirList(dirFlagArray & nameFlagArray);

% iterate over all subjects
for iSub = 1:numel(subDirList)
    
    subjectDir      = subDirList(iSub).name;
    sesDirList      = dir(fullfile(tempDir, subjectDir));
    
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
            sesEEGDir       = fullfile(tempDir, subjectDir, sesDir.name, 'eeg'); 
            sesBEHDir       = fullfile(tempDir, subjectDir, sesDir.name, 'beh'); 
            sesFilesEEG     = dir(sesEEGDir)';
            sesFilesBEH     = dir(sesBEHDir)';
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
        bidsModality    = bidsNameSplit{end}(1:end-4);                      % this string includes modality and extension
        extension       = bidsNameSplit{end}(end-3:end);
        
        if find(strncmp(bidsNameSplit,'ses',3))
            isMultiSession      = 1; 
            sessionName         = bidsNameSplit{strncmp(bidsNameSplit,'ses',3)}(5:end);
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
                bemobilName     = [bemobil_config.filename_prefix, num2str(subjectNr), '_' sessionName '_' bemobilModality, '_rec', runIndex, '_old', extension];
            else
                bemobilName     = [bemobil_config.filename_prefix, num2str(subjectNr), '_' sessionName '_' bemobilModality, '_old', extension];
            end
            
            data    = pop_loadset('filepath', dataDir, 'filename', bidsName);
            pop_saveset(data, 'filepath', newDir, 'filename', bemobilName);
            
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
nameArray       = nameArray(nameFlagArray); 
subDirList      = subDirList(nameFlagArray);

% select participant data to process
subFlagArray    = true(size(nameArray)); 
if ~isempty(numericalIDs)
    for iN = 1:numel(nameArray)
        subjectNr = str2double(nameArray{iN}(numel(bemobil_config.filename_prefix) + 1:end));
        if ~ismember(subjectNr,numericalIDs)
            subFlagArray(iN) = false; 
        end
    end
    subDirList = subDirList(subFlagArray); 
end

% merge all the run files if needed
for iSub = 1:numel(subDirList)
    
    % list all files in the subject folder
    subjectFiles = dir(fullfile(targetDir, subDirList(iSub).name));
    
    % iterate over sessions
    for iSes = 1:numel(bemobil_config.filenames)
        
        % check if final session file names are there - if so, skip
        subjectNr = str2double(subDirList(iSub).name(numel(bemobil_config.filename_prefix) + 1:end));
        sessionNameEEG = [bemobil_config.filename_prefix, num2str(subjectNr), '_', bemobil_config.filenames{iSes} '_EEG.set'];
       
        if isfile(fullfile(targetDir, subDirList(iSub).name, sessionNameEEG))
            disp(['File ' sessionNameEEG ' already found - skipping import for this session'])
            continue; 
        end
        
        % find all EEG data
        eegFiles = {subjectFiles(contains({subjectFiles.name}, [bemobil_config.filenames{iSes} '_EEG']) & contains({subjectFiles.name}, '_old.set')).name};
        eegFiles = natsortfiles(eegFiles); % not using natsortorder here - potentially problematic for more than 10 runs? (implausible)  
        
        if numel(eegFiles) > 1
            % if there are multiple runs, merge them all
            EEGMerged = pop_loadset('filepath',fullfile(targetDir, subDirList(iSub).name),'filename', eegFiles{1});
            for iFile = 2:numel(eegFiles)
                disp('merging files')
                [EEG2]      = pop_loadset('filepath', fullfile(targetDir, subDirList(iSub).name),'filename' ,eegFiles{iFile});
                EEGMerged   = pop_mergeset(EEGMerged, EEG2);
            end
            EEG                 = EEGMerged;
            EEGFileNameWithRun  = eegFiles{iFile}; 
            nameSplit           = regexp(EEGFileNameWithRun,'_', 'split'); % remove _rec entity
            nameJoined          = join(nameSplit(1:end-2),'_');
            EEGSessionFileName  = [nameJoined{1} '.set'];
        elseif numel(eegFiles) == 1
            EEG                 = pop_loadset('filepath',fullfile(targetDir, subDirList(iSub).name),'filename', eegFiles{1});
            EEGSessionFileName  = eegFiles{1};
        else
            warning(['No EEG file found in subject dir ' subDirList(iSub).name ', session ' bemobil_config.filenames{iSes}] )
        end
        
        % resample EEG
        %------------------------------------------------------------------
        if isempty(bemobil_config.resample_freq)
            newSRate        = round(EEG.srate);    
            warning(['No resample frequency specified - data is still resampled to the nearest integer ' num2str(newSRate) 'Hz'])
        else
            newSRate        = bemobil_config.resample_freq;
        end
        
        % Note that in fieldtrip time is in seconds
        newTimes                = (EEG.times(1):1000/newSRate:EEG.times(end))/1000;
        resamplecfg.time        = {newTimes};
        resamplecfg.detrend     = 'no';
        resamplecfg.extrapval   = nan;
        EEG.group = 1; EEG.condition = 1;
        ftData                  = eeglab2fieldtrip( EEG, 'raw', 'none' );
        resampledData           = ft_resampledata(resamplecfg, ftData);
        EEG.data                = resampledData.trial{1};
        EEG.srate               = newSRate;
        EEG.times               = newTimes*1000; % convert back to miliseconds
        EEG.pnts                = size(EEG.data,2);
        EEG.urevent             = EEG.event;
        for iE = 1:numel(EEG.event)
            EEG.event(iE).latency        = find(EEG.times > EEG.urevent(iE).latency,1,'first');
        end
        
        EEG.setname = EEG.filename(1:end-8); 
        
        % checkset 
        EEG = eeg_checkset(EEG); 
        
        % save merged EEG file for the session
        EEG = pop_saveset(EEG, 'filename',[EEGSessionFileName(1:end-8) EEGSessionFileName(end-3:end)],'filepath',fullfile(targetDir, subDirList(iSub).name));
        disp(['Saved session file ' EEGSessionFileName(1:end-8) EEGSessionFileName(end-3:end)])
        
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
            dataFiles = {subjectFiles(contains({subjectFiles.name}, [bemobil_config.filenames{iSes} '_' bemobilModality]) & contains({subjectFiles.name}, '_old.set')).name};
            
            if numel(dataFiles) > 1
                % if there are multiple runs, merge them all 
                DATAMerged =  pop_loadset('filepath', fullfile(targetDir, subDirList(iSub).name),'filename', dataFiles{1});
                for iFile = 2:numel(dataFiles)
                    DATA2           = pop_loadset('filepath', fullfile(targetDir, subDirList(iSub).name), 'filename' ,dataFiles{iFile});
                    DATAMerged      = pop_mergeset(DATAMerged, DATA2);
                end
                DATA                    = DATAMerged;
                DATAFileNameWithRun     = dataFiles{iFile};
                nameSplit               = regexp(DATAFileNameWithRun,'_', 'split'); % remove _run entity
                nameJoined              = join(nameSplit(1:end-2),'_');
                DATASessionFileName      = [nameJoined{1} '.set'];
            elseif numel(dataFiles) == 1
                DATA             = pop_loadset('filepath', fullfile(targetDir, subDirList(iSub).name), 'filename', dataFiles{1});
                DATASessionFileName  = dataFiles{1};
            else
                warning(['No file of modality ' bemobilModality ' found in subject dir ' subDirList(iSub).name ', session ' bemobil_config.filenames{iSes}] )
            end 
        end
       
        % resample DATA
        %------------------------------------------------------------------
        if isempty(bemobil_config.resample_freq)
            newSRate        = EEG.srate;
            warning(['No resample frequency specified -' bemobilModality 'data is still resampled to match EEG srate, ' num2str(newSRate) 'Hz'])
        else
            newSRate        = bemobil_config.resample_freq;
        end
        
        % unwrap any kind of angular data before resampling 
        angleind = []; 
        for Ci = 1:numel(DATA.chanlocs)
            if  contains(DATA.chanlocs(Ci).units, 'rad') 
                angleind =  [angleind Ci]; 
            end
        end
        
        % unwrap any kind of angular data before resampling 
        DATA.data(angleind,:)   = unwrap(DATA.data(angleind,:));
        
        % resample to EEG times 
        newTimes                = (EEG.times(1):1000/newSRate:EEG.times(end))/1000;
        resamplecfg.time        = {newTimes};
        resamplecfg.detrend     = 'no';
        resamplecfg.extrapval   = nan; 
        DATA.group = 1; DATA.condition = 1;
        ftData                  = eeglab2fieldtrip(DATA, 'raw', 'none' ); % here offset information is lost (onset corrected to zero): start time to be accounted for 
        ftData.time{1}          = ftData.time{1} + DATA.etc.starttime; % offset has to be added here 
        resampledData           = ft_resampledata(resamplecfg, ftData); % the resamplecfg used for EEG
        
        % assign the data back to the .set file 
        DATA.data                = resampledData.trial{1};
        
        % wrap back angular data to [pi, -pi]
        DATA.data(angleind,:)    = wrapToPi(DATA.data(angleind,:));
        DATA.srate               = newSRate;
        DATA.times               = resampledData.time{1}*1000;
        DATA.pnts                = size(DATA.data,2);
        
        % copy events from EEG
        DATA.event = EEG.event;
        
        % checkset 
        DATA = eeg_checkset(DATA);
        
        % save merged data file for the session
        DATA = pop_saveset(DATA, 'filename', [DATASessionFileName(1:end-8) DATASessionFileName(end-3:end)],'filepath',fullfile(targetDir, subDirList(iSub).name));
        disp(['Saved session file ' [DATASessionFileName(1:end-8) DATASessionFileName(end-3:end)]])
        
    end
    
    % remove unnecessary files prior to merging
    subjectFiles = dir(fullfile(targetDir, subDirList(iSub).name));
    toDelete = {subjectFiles(contains({subjectFiles.name}, '_old.fdt') | contains({subjectFiles.name}, '_old.set')).name};
    for iD = 1:numel(toDelete)
        delete(fullfile(targetDir, subDirList(iSub).name, toDelete{iD}));
    end
    
    if mergeAllSessions
        
%         % merge all EEG data over sessions 
%         for iFile = 1:numel(EEGFiles)
%             pop_mergeset
%         end
%         
%         % merge all other data over sessions
%         for iType = 1:numel(dataTypes)
%             for iFile = 1:numel(DataFiles)
%                 pop_mergeset
%             end
%         end
      
    end
    
end


end
