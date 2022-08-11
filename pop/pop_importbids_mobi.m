% pop_importbids_mobi() - Import BIDS format folder structure into an EEGLAB
%                    study.
% Usage:
%   >> [STUDY ALLEEG] = pop_importbids(bidsfolder);
%   >> [STUDY ALLEEG] = pop_importbids(bidsfolder, 'key', value);
%
% Inputs:
%   bidsfolder - a loaded epoched EEG dataset structure.
%     options are 'bidsevent', 'bidschanloc' of be turned 'on' (default) or 'off'
%                 'outputdir' default is bidsfolder/derivatives
%                 'studyName' default is eeg
%
% Optional inputs:
%  'studyName'   - [string] name of the STUDY
%  'bidsevent'   - ['on'|'off'] import events from BIDS .tsv file and
%                  ignore events in raw binary EEG files.
%  'bidschanloc' - ['on'|'off'] import channel location from BIDS .tsv file
%                  and ignore locations (if any) in raw binary EEG files.
%  'outputdir'   - [string] output folder (default is to use the BIDS
%                  folders).
%  'eventtype'   - [string] BIDS event column to use for EEGLAB event types.
%                  common choices are usually 'trial_type' or 'value'.
%                  Default is 'value'.
%  'datatypes'   - [cell array of strings] types of other data modalities than EEG
%                  current implementation only for motion according to BEP 029
%  'matchchanlocs' - [2x NChan array of strings] in case electrode names 
%                  in electrodes.tsv and channels.tsv do not match with
%                  each other. First column contains labels in electrodes.tsv
%                  and the second column contains lables in channels.tsv. 
%                  Resulting chanloc will take labels from eloc file. 
%                      example : {'g1', 'G01'; 'g2', 'G02'; ...}
%                  Use empty string for a missing chanloc 
%                      example : {'', 'N01'; 'n2', 'N02'; ...}
%  
% Outputs:
%   STUDY   - EEGLAB STUDY structure
%   ALLEEG  - EEGLAB ALLEEG structure
%   bids    - bids structure
%
% Authors: Arnaud Delorme, SCCN, INC, UCSD, January, 2019
%         Cyril Pernet, University of Edinburgh
%
% Example:
% pop_importbids('/data/matlab/bids_matlab/rishikesh_study/BIDS_EEG_meditation_experiment');

% Copyright (C) Arnaud Delorme, 2018
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 2 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program; if not, to the Free Software
% Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

function [STUDY, ALLEEG, bids, stats, commands] = pop_importbids_mobi(bidsFolder, varargin)

STUDY = [];
ALLEEG = [];
bids = [];
stats = [];
commands = '';
if nargin < 1
    bidsFolder = uigetdir('Pick a BIDS folder');
    if isequal(bidsFolder,0), return; end
    
    cb_select = [ 'tmpfolder = uigetdir;' ...
        'if ~isequal(tmpfolder, 0)' ...
        '   set(findobj(gcbf, ''tag'', ''folder''), ''string'', tmpfolder);' ...
        'end;' ...
        'clear tmpfolder;' ];
    type_fields = { 'value' 'trial_type' };
    
    cb_event = 'set(findobj(gcbf, ''userdata'', ''bidstype''), ''enable'', fastif(get(gcbo, ''value''), ''on'', ''off''));';
    promptstr    = { ...
        { 'style'  'text'       'string' 'Enter study name (default is BIDS folder name)' } ...
        { 'style'  'edit'       'string' '' 'tag' 'studyName' } ...
        {} ...
        { 'style'  'checkbox'   'string' 'Overwrite channel locations with BIDS channel location files' 'tag' 'chanlocs' 'value' 0 } ...
        { 'style'  'checkbox'  'string' 'Overwrite events with BIDS event files and use this BIDS field for event type' 'tag' 'events' 'value' 0 'callback' cb_event } ...
        { 'style'  'popupmenu'  'string' type_fields 'tag' 'typefield' 'value' 1 'userdata' 'bidstype'  'enable' 'off' } ...
        {} ...
        { 'style'  'text'       'string' 'Study output folder' } ...
        { 'style'  'edit'       'string' fullfile(bidsFolder, 'derivatives') 'tag' 'folder' 'HorizontalAlignment' 'left' } ...
        { 'style'  'pushbutton' 'string' '...' 'callback' cb_select } ...
        };
    geometry = {[2 1.5], 1, 1,[1 0.25],1,[1 2 0.5]};
    
    [~,~,~,res] = inputgui( 'geometry', geometry, 'geomvert', [1 0.5, 1 1 0.5 1], 'uilist', promptstr, 'helpcom', 'pophelp(''pop_importbids'')', 'title', 'Import BIDS data -- pop_importbids()');
    if isempty(res), return; end
    
    options = { 'eventtype' type_fields{res.typefield} };
    if res.events,    options = { options{:} 'bidsevent' 'on' };   else options = { options{:} 'bidsevent' 'off' }; end
    if res.chanlocs,  options = { options{:} 'bidschanloc' 'on' }; else options = { options{:} 'bidschanloc' 'off' }; end
    if ~isempty(res.folder),  options = { options{:} 'outputdir' res.folder }; end
    if ~isempty(res.studyName),  options = { options{:} 'studyName' res.studyName }; end
else
    options = varargin;
end

[~,defaultStudyName] = fileparts(bidsFolder);
opt = finputcheck(options, { ...
    'bidsevent'      'string'    { 'on' 'off' }    'on';  ...
    'bidschanloc'    'string'    { 'on' 'off' }    'on'; ...
    'metadata'       'string'    { 'on' 'off' }    'off'; ...
    'eventtype'      'string'    {  }              'value'; ...
    'outputdir'      'string'    { }               fullfile(bidsFolder,'derivatives'); ...
    'studyName'      'string'    { }               defaultStudyName; ...
    'datatypes'      'cell'      { }               {}; ...
    'matchchanlocs'  'cell'      { }               {}; ...
    'participants'   'integer'   []                [] ...  
    }, 'pop_importbids');
if isstr(opt), error(opt); end

% Options:
% - copy folder
% - use channel location and event

% load change file
changesFile = fullfile(bidsFolder, 'CHANGES');
bids.CHANGES = '';
if exist(changesFile,'File')
    bids.CHANGES = importalltxt( changesFile );
end

% load Readme file
readmeFile = fullfile(bidsFolder, 'README');
bids.README = '';
if exist(readmeFile,'File')
    bids.README = importalltxt( readmeFile );
end

% load dataset description file
dataset_descriptionFile = fullfile(bidsFolder, 'dataset_description.json');
bids.dataset_description = '';
if exist(dataset_descriptionFile,'File')
    bids.dataset_description = jsondecode(importalltxt( dataset_descriptionFile ));
end

% load participant tsv file
participantsFile = fullfile(bidsFolder, 'participants.tsv'); 
bids.participants = '';
if exist(participantsFile,'File')
    bids.participants = importtsv( participantsFile );
end

% load participant json file
participantsJSONFile = fullfile(bidsFolder, 'participants.json');
bids.participantsJSON = '';
if exist(participantsJSONFile,'File')
    bids.participantsJSON = jsondecode(importalltxt( participantsJSONFile ));
end

if isempty(bids.participants)
    participantFolders = dir(fullfile(bidsFolder, 'sub-*'));
    bids.participants = { participantFolders.name }';
end

% scan participants
count = 1;
commands = {};
task = [ 'task-' bidsFolder ];
bids.data = [];
inconsistentChannels = 0;
inconsistentEvents   = 0;

% check for data types 
datatypes = opt.datatypes;

for iSubject = 2:size(bids.participants,1)
    
    % check if participants list is given; if not, process all
    if ~isempty(opt.participants)
        
        % extract numerical participant ID 
        pNumeric = str2double(bids.participants{iSubject,1}(5:end)); 
        
        % if participant list is given, skip participants who are not listed
        if  ~ismember(pNumeric, opt.participants)
            continue;
        end
    end
    
    parentSubjectFolder = fullfile(bidsFolder   , bids.participants{iSubject,1});
    outputSubjectFolder = fullfile(opt.outputdir, bids.participants{iSubject,1});
    
    nSessions = 1; 
   
    % find folder containing eeg & load scans.tsv 
    if exist(fullfile(parentSubjectFolder, 'eeg'),'dir')
        isMultiSession = false; 
        subjectFolder = { fullfile(parentSubjectFolder, 'eeg') };
        subjectFolderOut = { fullfile(outputSubjectFolder, 'eeg') };
    else
        isMultiSession = true; 
        subFolders = dir(fullfile(parentSubjectFolder, 'ses*'));
        nSessions = length(subFolders);
        subjectFolder    = {};
        subjectFolderOut = {};
        
        for iFold = 1:nSessions
            subjectFolder{   iFold} = fullfile(parentSubjectFolder, subFolders(iFold).name, 'eeg');
            subjectFolderOut{iFold} = fullfile(outputSubjectFolder, subFolders(iFold).name, 'eeg');
            if ~exist(subjectFolder{iFold},'dir')
                subjectFolder{   iFold} = fullfile(parentSubjectFolder, subFolders(iFold).name, 'meg');
                subjectFolderOut{iFold} = fullfile(outputSubjectFolder, subFolders(iFold).name, 'meg');
            end
        end
       
    end
    
    % find folders containing data of other modalities than EEG
    for iType = 1:numel(datatypes)
        dataType            = datatypes{iType};
        
        if strcmpi(dataType, 'motion')
            modalityFolders     = [{'eeg'}, {'beh'}, {'motion'}]; % this may be dependent on the modality. For now we assume that the data is either in beh, motion or EEG     
        else 
            modalityFolders     = [{'eeg'}, {'beh'}]; 
        end 
        
        % iterate over modality specific folders to find the data files
        for iMod     =  1:numel(modalityFolders)
            if ~isMultiSession
                modFolderFullfile = fullfile(parentSubjectFolder, modalityFolders{iMod}); 
                if exist(modFolderFullfile,'dir')
                    filesInFolder       = dir(modFolderFullfile);
                    fileNamesInFolder   = {filesInFolder.name};
                    fileNamesInFolder(ismember(fileNamesInFolder,{'.','..'})) = [];
                    if any(strcmp(dataType, cellfun(@(x)x(end-length(dataType)-3:end-4), fileNamesInFolder, 'uniformoutput', false)))
                        disp(['Found ' dataType 'data in folder '  modalityFolders{iMod}])
                        % if the folder in which you found the data is not included in the search/output folder list, add them
                        if ~strcmp(modFolderFullfile, cellfun(@(x)x, subjectFolder, 'uniformoutput', false))
                            subjectFolder{end+1}    = modFolderFullfile; 
                            subjectFolderOut{end+1} = fullfile(outputSubjectFolder, modalityFolders{iMod});
                        end
                    end
                end
            else
                %subFolders = dir(fullfile(parentSubjectFolder, 'ses*'));
                for iFold = 1:nSessions
                    modFolderFullfile = fullfile(parentSubjectFolder, subFolders(iFold).name, modalityFolders{iMod});
                    if exist(modFolderFullfile,'dir')
                        filesInFolder       = dir(modFolderFullfile);
                        fileNamesInFolder   = {filesInFolder.name};
                        fileNamesInFolder(ismember(fileNamesInFolder,{'.','..'})) = [];
                        if any(strcmp(dataType, cellfun(@(x)x(end-length(dataType)-3:end-4), fileNamesInFolder, 'uniformoutput', false)))
                            disp(['Found ' dataType 'data in folder '  modalityFolders{iMod}])
                            % if the folder in which you found the data is not included in the search/output folder list, add them
                            if ~strcmp(modFolderFullfile, cellfun(@(x)x, subjectFolder, 'uniformoutput', false))
                                subjectFolder{end+1}    = modFolderFullfile; 
                                subjectFolderOut{end+1} = fullfile(outputSubjectFolder, subFolders(iFold).name, modalityFolders{iMod});
                            end
                        end
                    end
                end
            end
        end   
    end
    
    % import data
    for iFold = 1:length(subjectFolder) % scan sessions
        if ~exist(subjectFolder{iFold},'dir')
            fprintf(2, 'No data folder for subject %s session %s\n', bids.participants{iSubject,1}, subFolders(iFold).name);
        else
            
            % scans.tsv for time synch information
            %-------------------------------------
            try
                try
                    scansFile     = searchparent(fileparts(subjectFolder{iFold}), '*_scans.tsv');
                catch
                    % for some reason parent search not working - a quick workaround
                    scansFile       = searchparent(fileparts(fileparts(subjectFolder{iFold})), '*_scans.tsv');
                end
            catch
            end
            
            if exist('scansFile', 'var')
                useScans = true;
                scans = loadfile( scansFile.name, scansFile);
                bids.data = setallfields(bids.data, [iSubject-1,iFold,1], struct('scans', {scans}));
            else
                useScans = false;
            end
            
            if contains(subjectFolder{iFold}, 'eeg') || contains(subjectFolder{iFold}, 'meg')
                
                % which raw data - with folder inheritance
                eegFile     = searchparent(subjectFolder{iFold}, '*eeg.*');
                if isempty(eegFile)
                    eegFile     = searchparent(subjectFolder{iFold}, '*_meg.*');
                end
                
                infoFile      = searchparent(subjectFolder{iFold}, '*_eeg.json');
                channelFile   = searchparent(subjectFolder{iFold}, '*_channels.tsv');
                elecFile      = searchparent(subjectFolder{iFold}, '*_electrodes.tsv');
                eventFile     = searchparent(subjectFolder{iFold}, '*_events.tsv');
                
                try
                    eventDescFile = searchparent(subjectFolder{iFold}, '*_events.json');
                catch
                    eventDescFile = searchparent(bidsFolder, '*_events.json');
                end
                
                % raw EEG data
                allFiles = { eegFile.name };
                ind = strmatch( 'json', cellfun(@(x)x(end-3:end), allFiles, 'uniformoutput', false) );
                if ~isempty(ind)
                    eegFileJSON = allFiles(ind);
                    allFiles(ind) = [];
                end
                ind = strmatch( '.set', cellfun(@(x)x(end-3:end), allFiles, 'uniformoutput', false) ); % EEGLAB
                if ~isempty(ind)
                    eegFileRawAll  = allFiles(ind);
                elseif length(allFiles) == 1
                    eegFileRawAll  = allFiles;
                else
                    ind = strmatch( '.eeg', cellfun(@(x)x(end-3:end), allFiles, 'uniformoutput', false) ); % BVA
                    if isempty(ind)
                        ind = strmatch( '.edf', cellfun(@(x)x(end-3:end), allFiles, 'uniformoutput', false) ); % EDF
                        if isempty(ind)
                            ind = strmatch( '.bdf', cellfun(@(x)x(end-3:end), allFiles, 'uniformoutput', false) ); % BDF
                            if isempty(ind)
                                ind = strmatch( '.fif', cellfun(@(x)x(end-3:end), allFiles, 'uniformoutput', false) ); % FIF
                                if isempty(ind)
                                    ind = strmatch( '.gz', cellfun(@(x)x(end-2:end), allFiles, 'uniformoutput', false) ); % FIF
                                    if isempty(ind)
                                        fprintf(2, 'No EEG file found for subject %s\n', bids.participants{iSubject,1});
                                    end
                                end
                            end
                        end
                    end
                    eegFileRawAll  = allFiles(ind);
                end
            end
            
            allOtherDataFiles   = {}; 
            
            % find files of other datatypes
            for iType = 1:numel(datatypes)
                try
                    otherDataFile       = searchparent(subjectFolder{iFold}, ['*_' datatypes{iType} '.tsv']);
                    allOtherDataFiles   = [allOtherDataFiles { otherDataFile.name }];
                    if isempty(allOtherDataFiles)
                        disp(['No data of type ' datatypes{iType} ' found in ' subjectFolder{iFold}])
                    end
                catch
                    disp(['No data of type ' datatypes{iType} ' found in ' subjectFolder{iFold}])
                end
            end

            % skip most import if set file with no need for modification
            for iFile = 1:length(eegFileRawAll)

                eegFileRaw = eegFileRawAll{iFile};
                [~,tmpFileName,fileExt] = fileparts(eegFileRaw);
                eegFileRaw     = fullfile(subjectFolder{   iFold}, eegFileRaw);
                eegFileNameOut = fullfile(subjectFolderOut{iFold}, [ tmpFileName '.set' ]);
                
                % check for existence of the file before proceding
                if ~exist(eegFileRaw, 'file')
                    continue
                end
                
                % what is the run
                iRun = 1;
                ind = strfind(eegFileRaw, '_run-');
                if ~isempty(ind)
                    tmpEegFileRaw = eegFileRaw(ind(1)+5:end);
                    indUnder = find(tmpEegFileRaw == '_');
                    iRun = str2double(tmpEegFileRaw(1:indUnder(1)-1));
                    if isnan(iRun) || iRun == 0
                        iRun = str2double(tmpEegFileRaw(1:indUnder(1)-2)); % rare case run 5H in ds003190/sub-01/ses-01/eeg/sub-01_ses-01_task-ctos_run-5H_eeg.eeg
                        if isnan(iRun) || iRun == 0
                            error('Problem converting run information'); 
                        end
                    end
                end
                
                % JSON information file
                infoData = loadfile([ eegFileRaw(1:end-8) '_eeg.json' ], infoFile);
                bids.data = setallfields(bids.data, [iSubject-1,iFold,iFile], infoData);
                
                % optional coordinate system file
                try
                    coordFile = searchparent(subjectFolder{iFold}, '*_coordsystem.json');
                    coordData = loadfile(coordFile.name, coordFile);
                    bids.data = setallfields(bids.data, [iSubject-1,iFold,iFile], struct('coordinfo', { coordData }));
                catch
                    disp('No coordinate system file found')
                end
                
                % extract task name
                underScores = find(tmpFileName == '_');
                if ~strcmpi(tmpFileName(underScores(end)+1:end), 'eeg')
                    if ~strcmpi(tmpFileName(underScores(end)+1:end), 'meg.fif')
                        if ~strcmpi(tmpFileName(underScores(end)+1:end), 'meg')
                            error('Data file name does not contain eeg or meg'); % theoretically impossible
                        end
                    end
                end
                
                if isempty(findstr('ses', tmpFileName(underScores(end-1)+1:underScores(end)-1)))
                    task = tmpFileName(underScores(end-1)+1:underScores(end)-1);
                end
                
                if ~strcmpi(fileExt, '.set') || strcmpi(opt.bidsevent, 'on') || strcmpi(opt.bidschanloc, 'on') || ~strcmpi(opt.outputdir, bidsFolder)
                    fprintf('Importing file: %s\n', eegFileRaw);
                    switch lower(fileExt)
                        case '.set' % do nothing
                            if strcmpi(opt.metadata, 'on')
                                EEG = pop_loadset( 'filename', eegFileRaw, 'loadmode', 'info' );
                            else
                                EEG = pop_loadset( 'filename', eegFileRaw );
                            end
                        case {'.bdf','.edf'}
                            EEG = pop_biosig( eegFileRaw ); % no way to read meta data only (because events in channel)
                        case '.eeg'
                            [tmpPath,tmpFileName,~] = fileparts(eegFileRaw);
                            if exist(fullfile(tmpPath, [tmpFileName '.vhdr']), 'file'), ext = '.vhdr'; else ext = '.VMRK'; end
                            if strcmpi(opt.metadata, 'on')
                                EEG = pop_loadbv( tmpPath, [tmpFileName ext], [], [], true );
                            else
                                EEG = pop_loadbv( tmpPath, [tmpFileName ext] );
                            end
                        case '.fif'
                            EEG = pop_fileio(eegFileRaw); % fif folder
                        case '.gz'
                            gunzip(eegFileRaw);
                            EEG = pop_fileio(eegFileRaw(1:end-3)); % fif folder
                        case '.ds'
                            EEG = pop_fileio(eegFileRaw); % fif folder
                        otherwise
                            error('No EEG data found for subject/session %s', subjectFolder{iFold});
                    end
                    EEGnodata = EEG;
                    EEGnodata.data = [];
                    bids.data = setallfields(bids.data, [iSubject-1,iFold,iFile], struct('EEG', EEGnodata));
               
                    % channel location data
                    % ---------------------
                    channelData = loadfile([ eegFileRaw(1:end-8) '_channels.tsv' ], channelFile);
                    elecData    = loadfile([ eegFileRaw(1:end-8) '_electrodes.tsv' ], elecFile);
                    bids.data = setallfields(bids.data, [iSubject-1,iFold,iFile], struct('chaninfo', { channelData }));
                    bids.data = setallfields(bids.data, [iSubject-1,iFold,iFile], struct('elecinfo', { elecData }));
                    if strcmpi(opt.bidschanloc, 'on')
                        
                        if isfield(bids.data(iSubject-1,iFold,iFile),'EEGReference') && ~isempty(bids.data(iSubject-1,iFold,iFile).EEGReference)
                            specified_ref = bids.data(iSubject-1,iFold,iFile).EEGReference;
                        else
                            specified_ref = 'n/a';
                        end
                        
                        % scan elecData and set aside other elocs that are
                        % not present in EEG data 
                        if ~isempty(elecData)
                            nameCol         = strcmp(elecData(1,:), 'name');
                            otherLocInd      = [];
                            for iE = 2:size(elecData,1)
                                switch  cell2mat(lower(elecData(iE,nameCol)))
                                    case lower(specified_ref)
                                        otherLocInd = [otherLocInd iE];
                                    case 'ref'
                                        otherLocInd = [otherLocInd iE];
                                    case 'nasion'
                                        otherLocInd = [otherLocInd iE];
                                    case 'leftear'
                                        otherLocInd = [otherLocInd iE];
                                    case 'rightear'
                                        otherLocInd = [otherLocInd iE];
                                end
                            end
                            EEG.etc.extralocs         = elecData(otherLocInd,:);
                            elecData(otherLocInd, :)  = [];
                        end
                        
                        if ~isempty(opt.matchchanlocs)
                            for iChan = 2:size(channelData,1)
                                % the fields below are all required
                                chanlocs(iChan-1).labels = channelData{iChan,1};
                                chanlocs(iChan-1).type   = channelData{iChan,2};
                                chanlocs(iChan-1).unit   = channelData{iChan,3};
                                if size(channelData,2) > 3
                                    chanlocs(iChan-1).status = channelData{iChan,4};
                                end
                                
                                if ~isempty(elecData) 
                                    elecNameCol         = find(strcmp(elecData(1,:),'name')); 
                                    matchRow            = find(strcmpi(opt.matchchanlocs(:,2), chanlocs(iChan-1).labels));
                                    
                                    if isempty(opt.matchchanlocs{matchRow,1})
                                        chanlocs(iChan-1).labels    = ['E' num2str(iChan-1)];
                                        chanlocs(iChan-1).X         = [];
                                        chanlocs(iChan-1).Y         = [];
                                        chanlocs(iChan-1).Z         = [];
                                    else
                                        iElec               = find(strcmpi(elecData(:,elecNameCol), opt.matchchanlocs{matchRow,1}));
                                        chanlocs(iChan-1).labels    = elecData{iElec,1};
                                        chanlocs(iChan-1).X         = elecData{iElec,2};
                                        chanlocs(iChan-1).Y         = elecData{iElec,3};
                                        chanlocs(iChan-1).Z         = elecData{iElec,4};
                                    end
                                    
                                end
                            end
                        else
                            chanlocs = [];
                            for iChan = 2:size(channelData,1)
                                % the fields below are all required
                                chanlocs(iChan-1).labels = channelData{iChan,1};
                                chanlocs(iChan-1).type   = channelData{iChan,2};
                                chanlocs(iChan-1).unit   = channelData{iChan,3};
                                if size(channelData,2) > 3
                                    chanlocs(iChan-1).status = channelData{iChan,4};
                                end
                                
                                if ~isempty(elecData) && iChan <= size(elecData,1)
                                    chanlocs(iChan-1).labels    = elecData{iChan,1};
                                    chanlocs(iChan-1).X         = elecData{iChan,2};
                                    chanlocs(iChan-1).Y         = elecData{iChan,3};
                                    chanlocs(iChan-1).Z         = elecData{iChan,4};
                                end
                            end
                        end
                        if length(chanlocs) ~= EEG.nbchan
                            warning('Different number of channels in channel location file and EEG file');
                            % check if the difference is due to non EEG channels
                            % list here https://bids-specification.readthedocs.io/en/stable/04-modality-specific-files/03-electroencephalography.html
                            keep = {'EEG','EOG','HEOG','VEOG'}; % keep all eeg related channels
                            tsv_eegchannels  = arrayfun(@(x) sum(strcmpi(x.type,keep)),chanlocs,'UniformOutput',true);
                            tmpchanlocs = chanlocs; tmpchanlocs(tsv_eegchannels==0)=[]; % remove non eeg related channels
                            chanlocs = tmpchanlocs; clear tmpchanlocs
                            disp('saving non-EEG channels ')
                        end
                        
                        if length(chanlocs) ~= EEG.nbchan
                            error('channel location file and EEG file do not have the same number of channels');
                        end
                        
                        if isfield(chanlocs, 'X')
                            EEG.chanlocs = convertlocs(chanlocs, 'cart2all');
                        else
                            EEG.chanlocs = chanlocs;
                        end
                        
                        % optional coordinate system file
                        try
                            coordData = loadfile([], coordFile);
                            EEG.etc.coordsys = coordData; 
                        catch
                            disp('No coordinate system file found')
                        end
                        
                    end
                    
                    % event data
                    % ----------
                    eventData = loadfile( [ eegFileRaw(1:end-8) '_events.tsv' ], eventFile);
                    bids.data = setallfields(bids.data, [iSubject-1,iFold,iFile], struct('eventinfo', {eventData}));
                    eventDesc = loadfile( [ eegFileRaw(1:end-8) '_events.json' ], eventDescFile);
                    bids.data = setallfields(bids.data, [iSubject-1,iFold,iFile], struct('eventdesc', {eventDesc}));
                    if strcmpi(opt.bidsevent, 'on') && ~isempty(eventData)                      
                        events = struct([]);
                        indTrial = strmatch( opt.eventtype, lower(eventData(1,:)), 'exact');
                        for iEvent = 2:size(eventData,1)
                            events(end+1).latency  = eventData{iEvent,1}*EEG.srate+1; % convert to samples
                            events(end).duration   = eventData{iEvent,2}*EEG.srate;   % convert to samples
                            if ~isempty(indTrial)
                                events(end).type = eventData{iEvent,indTrial};
                            end
                            for iField = 1:length(eventData(1,:))
                                if ~strcmpi(eventData{1,iField}, 'onset') && ~strcmpi(eventData{1,iField}, 'duration')
                                    events(end).(eventData{1,iField}) = eventData{iEvent,iField};
                                end
                            end
                        end
                        EEG.event = events;
                        
                        % rename event fields to account for different usage of some field names 
                        if strcmpi(EEG.event(1).type, 'Markers') 
                           [EEG.event(:).type] = EEG.event(:).value;
                            EEG.event       = rmfield(EEG.event, 'value'); 
                            EEG.event       = rmfield(EEG.event, 'sample');
                        end
                        
                        EEG = eeg_checkset(EEG, 'eventconsistency');
                    end
                    
                    % copy information inside dataset
                    EEG.subject = bids.participants{iSubject,1};
                    EEG.session = iFold;
                    EEG.etc.nominal_srate = infoData.SamplingFrequency;
                    if strcmpi(opt.metadata, 'off')
                        if exist(subjectFolderOut{iFold},'dir') ~= 7
                            mkdir(subjectFolderOut{iFold});
                        end
                        EEG = pop_saveset( EEG, eegFileNameOut);
                    end
                end
       
                count = count+1;
                
                % check dataset consistency
                bData = bids.data(iSubject-1,iFold,iFile);
                if ~isempty(bData.chaninfo)
                    if size(bData.chaninfo,1)-1 ~= bData.EEG.nbchan
                        fprintf(2, 'Warning: inconsistency detected, %d channels in BIDS file vs %d in EEG file for %s\n', size(bData.chaninfo,1)-1, bData.EEG.nbchan, [tmpFileName,fileExt]);
                        inconsistentChannels = inconsistentChannels+1;
                    end
                end
                if ~isempty(bData.eventinfo)
                    if size(bData.eventinfo,1)-1 ~= length(bData.EEG.event)
                        fprintf(2, 'Warning: inconsistency detected, %d events in BIDS file vs %d in EEG file for %s\n', size(bData.eventinfo,1)-1, length(bData.EEG.event), [tmpFileName,fileExt]);
                        inconsistentEvents = inconsistentEvents+1;
                    end
                end  
                
                % check dataset consistency
                bData = bids.data(iSubject-1,iFold,iFile);
                if ~isempty(bData.chaninfo)
                    if size(bData.chaninfo,1)-1 ~= bData.EEG.nbchan
                        fprintf(2, 'Warning: inconsistency detected, %d channels in BIDS file vs %d in EEG file for %s\n', size(bData.chaninfo,1)-1, bData.EEG.nbchan, [tmpFileName,fileExt]);
                        inconsistentChannels = inconsistentChannels+1;
                    end
                end
                
                if ~isempty(bData.eventinfo)
                    if size(bData.eventinfo,1)-1 ~= length(bData.EEG.event)
                        fprintf(2, 'Warning: inconsistency detected, %d events in BIDS file vs %d in EEG file for %s\n', size(bData.eventinfo,1)-1, length(bData.EEG.event), [tmpFileName,fileExt]);
                        inconsistentEvents = inconsistentEvents+1;
                    end
                end
            end % end for eegFileRaw
           
            for iFile = 1:length(allOtherDataFiles)
       
                otherFileRaw    = allOtherDataFiles{iFile};
                bids.otherdata  = [];
                
                % extract data type
                splitName       = regexp(otherFileRaw,'_','split');
                datatype        = splitName{end}(1:end-4);
             
                joinedName          = join(splitName, '_'); 
                otherFileJSON       = [joinedName{1}(1:end-3) 'json']; 
                otherFileChannels   = [joinedName{1}(1:end-10) 'channels.tsv'];
                
                % check needed files according to the data type 
                switch datatype
                    case 'motion'
                        
                        % JSON information file
                        infoFileMotion  = searchparent(subjectFolder{iFold}, otherFileJSON);
                        infoData        = loadfile(infoFileMotion.name, infoFileMotion);
                        bids.otherdata  = setallfields(bids.otherdata, [iSubject-1,iFold,iFile], infoData);
                        
                        % channel file 
                        importChan          = true;
                        channelFileMotion   = searchparent(subjectFolder{iFold}, otherFileChannels);
                        channelData         = loadfile(channelFileMotion(1).name, channelFileMotion); % this might have to change along with motion BEP
                                      
                        % coordinate system file (potentially shared with EEG)
                        importCoord 	= false;
                        
                    case 'physio'
                        
                        % JSON information file
                        infoFilePhysio  = searchparent(subjectFolder{iFold}, otherFileJSON);
                        infoData        = loadfile(infoFilePhysio.name, infoFilePhysio);
                        bids.otherdata  = setallfields(bids.otherdata, [iSubject-1,iFold,iFile], infoData);
                        
                        % channel file (for physio data, hidden in json file)
                        importChan = true;
                        channelData = {'name', 'type', 'units'}; 
                        for Ci = 1:numel(infoData.Columns)
                            channelData{end + 1,1} = [char(infoData.Columns{Ci})] ;
                        end
                        
                        % coordinate file (none)
                        importCoord = false;
                        
                    otherwise
                        importChan = false;
                        disp(['No channel file for undefined data type ' datatype ])
                end
                
                % construct full file paths
                [~,tmpFileName,fileExt] = fileparts(otherFileRaw);
                otherFileRaw     = fullfile(subjectFolder{   iFold}, otherFileRaw);
                otherFileNameOut = fullfile(subjectFolderOut{iFold}, [ tmpFileName '.set' ]);
                 
                % check for existence of the file before proceding
                if ~exist(otherFileRaw, 'file')
                    continue
                end
                
                % what is the run
                iRun = 1;
                ind = strfind(otherFileRaw, '_run-');
                if ~isempty(ind)
                    tmpOtherFileRaw = otherFileRaw(ind(1)+5:end);
                    indUnder = find(tmpOtherFileRaw == '_');
                    iRun = str2double(tmpOtherFileRaw(1:indUnder(1)-1));
                    if isnan(iRun) || iRun == 0
                        iRun = str2double(tmpOtherFileRaw(1:indUnder(1)-2)); % rare case run 5H in ds003190/sub-01/ses-01/eeg/sub-01_ses-01_task-ctos_run-5H_eeg.eeg
                        if isnan(iRun) || iRun == 0
                            error('Problem converting run information');
                        end
                    end
                end
                
                
                if ~strcmpi(fileExt, '.set') || strcmpi(opt.bidsevent, 'on') || strcmpi(opt.bidschanloc, 'on') || ~strcmpi(opt.outputdir, bidsFolder)
                    fprintf('Importing file: %s\n', otherFileRaw);
                    switch lower(fileExt)
                        case '.set' % do nothing
                            if strcmpi(opt.metadata, 'on')
                                DATA = pop_loadset( 'filename', otherFileRaw, 'loadmode', 'info' );
                            else
                                DATA = pop_loadset( 'filename', otherFileRaw );
                            end
                        case '.tsv'
                            DATA        = eeg_emptyset;
                            
                            % check if the .tsv file has a line of header
                            fid = fopen(otherFileRaw);
                            header = textscan(fid, '%s',1, 'delimiter', '\n');
                            header = header{1};
                            headers = regexp(header{1},'\t','split');
                            fclose(fid);
                            
                            if isempty(header)
                                DATA.data   = dlmread(otherFileRaw,'\t')';
                            else
                                DATA.data   = dlmread(otherFileRaw,'\t',1,0)';
                            end
                            
                            if strcmp(datatype,'motion') 
                                DATA.srate                  = infoData.SamplingFrequencyEffective; 
                                DATA.etc.nominal_srate      = infoData.SamplingFrequency;
                            else
                                try
                                    DATA.srate  = infoData.SamplingFrequencyEffective; % Actual sampling rate used in motion data. Note that the unit of the time must be in second.
                                catch
                                    DATA.srate  = infoData.SamplingFrequency; % Generic physio data 
                                end
                            end
                            
                            useLatency = 0;
                            if strcmp(datatype,'motion')
                                latencyInd = find(strcmpi(channelData(:,strcmp(channelData(1,:),'type')), 'latency'));
                                useLatency = ~isempty(latencyInd);
                                if useLatency
                                    latencyHeader   = channelData{latencyInd,strcmp(channelData(1,:),'name')};
                                    latencyRowInData = find(strcmp(headers, latencyHeader));
                                end
                            elseif strcmp(datatype,'physio')
                                % check if the tracking system comes with latency
                                latencyInd      = find(contains(channelData(:,strcmp(channelData(1,:),'name')), 'latency'));
                                useLatency = ~isempty(latencyInd);
                                if useLatency
                                    latencyHeader   = channelData{latencyInd,strcmp(channelData(1,:),'name')};
                                    latencyRowInData = find(strcmp(headers, latencyHeader));
                                end
                            end
                            
                            % reconstruct time : use scans.tsv for synching
                            % it computes offset between motion and eeg data
                            if useScans
                                
                                scansData   = bids.data(iSubject-1,iFold,1).scans; 
                                for Coli = 1:size(scansData, 2)
                                    if strcmp(scansData{1,Coli}, 'acq_time') 
                                        acqTimeColi = Coli;  
                                    elseif strcmp(scansData{1,Coli}, 'filename')
                                        fNameColi = Coli; 
                                    end
                                end 
                                
                                for Rowi = 1:size(scansData, 1)
                                    
                                    sesString = ''; 
                                    taskString = ''; 
                                    runString = ''; 
                                    trackSysString = ''; 
                                    
                                    if exist('tracksys', 'var') && ~isempty(tracksys)
                                        trackSysString = tracksys; 
                                    end
                                        
                                    for SNi = 1:numel(splitName) 
                                        if contains(splitName{SNi}, 'ses-')
                                            sesString = splitName{SNi}(5:end);
                                        elseif contains(splitName{SNi}, 'task-')
                                            taskString = splitName{SNi}(6:end);
                                        elseif contains(splitName{SNi}, 'run-')
                                            runString = splitName{SNi}(5:end);
                                        end
                                    end
                                    
                                    [~,rawName] = fileparts(otherFileRaw);  
                                    % find files that matches in session, task, tracking system (in case it is motion data), and run
                                    if contains(scansData{Rowi,fNameColi}, 'eeg.') &&...
                                            contains(scansData(Rowi,fNameColi), sesString) && contains(scansData(Rowi,fNameColi), taskString) &&...
                                            contains(scansData(Rowi,fNameColi), runString)
                                        eegAcqTime      = scansData(Rowi,acqTimeColi); 
                                    elseif contains(scansData(Rowi,fNameColi), sesString) && contains(scansData(Rowi,fNameColi), taskString) &&...
                                           contains(scansData(Rowi,fNameColi), runString) && contains(scansData(Rowi,fNameColi), trackSysString) &&...
                                           contains(scansData(Rowi,fNameColi), datatype)
                                        otherAcqTime    = scansData(Rowi,acqTimeColi); 
                                    end
                                end
                                
                                startTime = seconds(datetime(otherAcqTime{1}, 'InputFormat', 'yyyy-MM-dd''T''HH:mm:ss.SSS') - datetime(eegAcqTime{1}, 'InputFormat', 'yyyy-MM-dd''T''HH:mm:ss.SSS'));
                                 
                            else
                                if isfield(infoData, 'StartTime')
                                    if isnumeric(infoData.StartTime)
                                        startTime  = infoData.StartTime;
                                    else
                                        startTime = 0;
                                        disp('Field Start time in motion json file is non-numeric - assume no offset to eeg data')
                                    end
                                else
                                    startTime = 0;
                                    disp('Field Start time in motion json file is empty - assume no offset to eeg data')
                                end
                            end
                            
                            DATA.etc.starttime = startTime;
                            
                            if useLatency
                                DATA.times = (DATA.data(latencyRowInData,:) - DATA.data(latencyRowInData,1))*1000;
                                DATA.data(latencyRowInData,:)   = [];
                            else
                                if isfield(infoData, 'TrackingSystems')
                                    % tsi should have been defined above when srate was being read in
                                    DATA.times  = (0:1000/infoData.TrackingSystems(tsi).SamplingFrequencyEffective:infoData.TrackingSystems(tsi).RecordingDuration*1000); % time is in ms
                                elseif isfield(infoData, 'RecordingDuration')
                                    DATA.times  = (0:1000/infoData.SamplingFrequency:infoData.RecordingDuration*1000); % time is in ms
                                else
                                    DATA.times  = (0:1000/infoData.SamplingFrequency:(size(DATA.data,2)/infoData.SamplingFrequency)*1000); % time is in ms
                                end
                            end
                            
                            DATA.nbchan = size(DATA.data,1);
                            DATA.pnts   = size(DATA.data,2);
                            
                        otherwise
                            error(['No ' datatype 'data found for subject/session ' subjectFolder{iFold}]);
                            
                    end
                    
                    bids.otherdata = setallfields(bids.otherdata, [iSubject-1,iFold,iFile], struct(datatype, []));                   
                    
                    if importChan
                        bids.otherdata = setallfields(bids.otherdata, [iSubject-1,iFold,iFile], struct('chaninfo', { channelData }));
                        if strcmpi(opt.bidschanloc, 'on')
                            chanlocs = [];
                            for iChan = 2:size(channelData,1)
                                % the fields below are all required
                                for iField = 1:size(channelData,2)
                                    fName = string(channelData(1,iField));
                                    if strcmp(fName, 'name')
                                        fName = 'labels';
                                    end
                                    chanlocs(iChan-1).(fName)   = channelData{iChan,iField};
                                end
                                if size(channelData,2) > 3
                                    chanlocs(iChan-1).status = channelData{iChan,4};
                                end
                            end
                            if useLatency
                                chanlocs(latencyRowInData)   = [];
                            end
                        end
                    end
                    
                    if importCoord
                        bids.otherdata = setallfields(bids.otherdata, [iSubject-1,iFold,iFile], struct('coordinfo', { coordData }));
                        DATA.etc.coordsys = coordData; 
                    end
                    
                    % copy information inside dataset
                    DATA.subject = bids.participants{iSubject,1};
                    DATA.session = iFold;
                    DATA.chanlocs = chanlocs; 
                    
                    if strcmpi(opt.metadata, 'off')
                        if exist(subjectFolderOut{iFold},'dir') ~= 7
                            mkdir(subjectFolderOut{iFold});
                        end
                        DATA = pop_saveset( DATA, otherFileNameOut);
                    end
                end
            end % end of other file Raw
        end
    end
end

% update statistics
% -----------------
% compute basic statistics
stats.README             = 0;
stats.TaskDescription    = 0;
stats.Instructions       = 0;
stats.EEGReference       = 0;
stats.PowerLineFrequency = 0;
stats.ChannelTypes       = 0;
stats.ElectrodePositions = 0;
stats.ParticipantsAgeAndSex = 0;
stats.SubjectArtefactDescription = 0;
stats.eventConsistency   = 0;
stats.channelConsistency = 0;
stats.EventDescription    = 0;
if ~isempty(bids.README), stats.README = 1; end
if isempty(bids.participants)
    error('No BIDS files found for specified participants')
end
if ismember('age'   , bids.participants(1,:)) && ismember('sex', bids.participants(1,:))
    stats.ParticipantsAgeAndSex = 1; 
end
if checkBIDSfield(bids, 'TaskDescription'),            stats.TaskDescription = 1; end
if checkBIDSfield(bids, 'Instructions'),               stats.Instructions = 1; end
if checkBIDSfield(bids, 'EEGReference'),               stats.EEGReference = 1; end
if checkBIDSfield(bids, 'PowerLineFrequency'),         stats.PowerLineFrequency = 1; end
if checkBIDSfield(bids, 'elecinfo'),                   stats.ElectrodePositions = 1; end
if checkBIDSfield(bids, 'eventdesc'),                  stats.EventDescription   = 1; end
if checkBIDSfield(bids, 'SubjectArtefactDescription'), stats.SubjectArtefactDescription   = 1; end
if isfield(bids.data, 'chaninfo') && ~isempty(bids.data(1).chaninfo) && ~isempty(strmatch('type', lower(bids.data(1).chaninfo(1,:)), 'exact'))
    stats.ChannelTypes = 1;
end
stats.channelConsistency = fastif(inconsistentChannels > 0, 0, 1);
stats.eventConsistency   = fastif(inconsistentEvents   > 0, 0, 1);



% check BIDS data field present
% -----------------------------
function res = checkBIDSfield(bids, fieldName)
res = false;
if isfield(bids.data, fieldName)
    fieldContent = { bids.data.(fieldName) };
    fieldContent(cellfun(@isempty, fieldContent)) = [];
    if ~isempty(fieldContent), res = true; end
end

% Import full text file
% ---------------------
function str = importalltxt(fileName)

str = [];
fid =fopen(fileName, 'r');
while ~feof(fid)
    str = [str 10 fgetl(fid) ];
end
str(1) = [];

% search parent folders (outward search) for the file of given fileName
% ---------------------
function outFile = searchparent(folder, fileName)
% search nestedly outward
% only get exact match and filter out hidden file
outFile = '';
parent = folder;
while ~any(arrayfun(@(x) strcmp(lower(x.name),'dataset_description.json'), dir(parent))) && isempty(outFile) % README indicates root BIDS folder
    outFile = filterHiddenFile(folder, dir(fullfile(parent, fileName)));
    parent = fileparts(parent);
end
if isempty(outFile)
    outFile = filterHiddenFile(folder, dir(fullfile(parent, fileName)));
end

function fileList = filterHiddenFile(folder, fileList)
isGoodFile = true(1,numel(fileList));
% loop to identify hidden files
for iFile = 1:numel(fileList) %'# loop only non-dirs
    % on OSX, hidden files start with a dot
    isGoodFile(iFile) = logical(~strcmp(fileList(iFile).name(1),'.'));
    if isGoodFile(iFile) && ispc
        % check for hidden Windows files - only works on Windows
        [~,stats] = fileattrib(fullfile(folder,fileList(iFile).name));
        if stats.hidden
            isGoodFile(iFile) = false;
        end
    end
end

% remove bad files
fileList = fileList(isGoodFile);

% import JSON or TSV file
function data = loadfile(localFile, globalFile)
[~,~,ext] = fileparts(localFile);
data = [];
localFile = dir(localFile);
if ~isempty(localFile)
    if strcmpi(ext, '.tsv')
        data = importtsv( fullfile(localFile(1).folder, localFile(1).name));
    else
        data = jsondecode( importalltxt( fullfile(localFile(1).folder, localFile(1).name) ));
    end        
elseif ~isempty(globalFile)
    if strcmpi(ext, '.tsv')
        data = importtsv( fullfile(globalFile(1).folder, globalFile(1).name));
    else
        data = jsondecode( importalltxt( fullfile(globalFile(1).folder, globalFile(1).name) ));
    end
end

% set structure
function sdata = setallfields(sdata, indices, newdata)
if isempty(newdata), return; end
if ~isstruct(newdata), error('Can only assign structures'); end
if length(indices) < 3, error('Must have 3 indices'); end
allFields = fieldnames(newdata);
for iField = 1:length(allFields)
    sdata(indices(1), indices(2), indices(3)).(allFields{iField}) = newdata.(allFields{iField});
end

% Import tsv file
% ---------------
function res = importtsv( fileName)

res = loadtxt( fileName, 'verbose', 'off', 'delim', 9);

for iCol = 1:size(res,2)
    % search for NaNs in numerical array
    indNaNs = cellfun(@(x)strcmpi('n/a', x), res(:,iCol));
    if ~isempty(indNaNs)
        allNonNaNVals = res(find(~indNaNs),iCol);
        allNonNaNVals(1) = []; % header
        testNumeric   = cellfun(@isnumeric, allNonNaNVals);
        if all(testNumeric)
            res(find(indNaNs),iCol) = { NaN };
        elseif ~all(~testNumeric)
            % Convert numerical value back to string
            res(:,iCol) = cellfun(@num2str, res(:,iCol), 'uniformoutput', false);
        end
    end
end
