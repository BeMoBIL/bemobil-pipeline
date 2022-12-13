function bemobil_bids2set(config)
% This function reads in BIDS datasets using the eeglab plugin
% "bids-matlab-tools" and reorganizes the output to be compatible with
% BeMoBIL pipeline.
%
% Usage
%       bemobil_bids2set(config);
%
% In
%       config.bids_target_folder     = 'P:\...SPOT_rotation\1_BIDS-data';  % required
%       config.set_folder             = 'P:\...SPOT_rotation\2_EEGlab-basic'; % required
%       config.subject                = 1;                                  % required, can also be an array in case importing a whole data set (e.g., config.subject = [1:10]; )
%       config.session_names          = {'body', 'joy'};                    % required, enter task name as a string, or enter a cell array when there are multiple sessions in the data set
%       config.overwrite              = 'on';                               % optional, default value 'off'
%       config.filename_prefix        = 'sub-';                             % optional, default value 'sub-'
%       config.resample_freq          = 500;                                % optional, default is 250 Hz
%       config.match_electrodes_channels = {'g1', 'G01'; 'g2', 'G02';...};  % optional, 2x NChan (number of channels in EEG data) array of strings, in case electrode names
%                                                                             in electrodes.tsv and channels.tsv do not match with
%                                                                             each other. First column contains labels in electrodes.tsv
%                                                                             and the second column contains lables in channels.tsv.
%                                                                             Resulting chanloc will take labels from eloc file.
%                                                                             Use empty string for a missing chanloc
%                                                                                   example : {'', 'N01'; 'n2', 'N02'; ...}
%       config.other_data_types        = {'motion','physio'};               % optional, default value is {}, required
%                                                                               if you want to load non-EEG data. Only
%                                                                               add the data type you have in your BIDS
%                                                                               data.
%       config.use_nominal_srate       = {'HTCVive','PhaseSpace'...};       % optional, in case you want to rely on nominal srate rather than time stamps for time synchronization, enter the names of corresponding tracksys names from motion bids files
%       
%
% Out
%       none
%       reorganizes data on disk
%
% required plugins
%       bva-io for brain vision data :
%               https://github.com/arnodelorme/bva-io
%
% author : seinjeung@gmail.com
%--------------------------------------------------------------------------

% input check and default value assignment
%--------------------------------------------------------------------------
% required fields
config = checkfield(config, 'bids_target_folder', 'required', '');
config = checkfield(config, 'set_folder', 'required', '');
config = checkfield(config, 'subject', 'required', '');
config = checkfield(config, 'session_names', 'required', '');

% default values
config = checkfield(config, 'filename_prefix', 'sub-', 'sub-');
config = checkfield(config, 'merged_filename', 'merged', 'merged');
config = checkfield(config, 'other_data_types', {}, '{}');
config = checkfield(config, 'resample_freq', 250, '250 Hz');
config = checkfield(config, 'overwrite', 'off', 'off');
config = checkfield(config, 'match_electrodes_channels', {}, 'none');
config = checkfield(config, 'use_nominal_srate', {}, 'none');
config = checkfield(config, 'do_aa_filter', 1, '1');

% check session input
if ~iscell(config.session_names)
    config.session_names = {config.session_names}; % convert string to cell in unisession case
end

% other data types than eeg (use BIDS suffix)
otherDataTypes  = config.other_data_types;

% construct the bids data directory
bidsDir         = fullfile(config.bids_target_folder);

% all runs and sessions are merged by default
targetDir       = fullfile(config.set_folder);
tempDir         = fullfile(targetDir, 'temp_bids');

if isfolder(tempDir)
    warning('Previous import seems to have been interrupted - removing temp_bids folder')
    rmdir(tempDir, 's')
end

numericalIDs        = config.subject;

if strcmp(config.overwrite, 'off')

    skipFlag            = false(size(numericalIDs));

    % check if final session file names are there - if so, skip
    for IDi = 1:numel(numericalIDs)
        [sesfilePath, sesfileName] = sessionfilename(targetDir, 'EEG', config, 1, numericalIDs(IDi));
        if isfile(fullfile(sesfilePath, sesfileName))
            disp(['File ' sesfileName ' found - skipping import for this participant'])
            skipFlag(IDi) = true;
        end
    end

    numericalIDs = numericalIDs(~skipFlag);

    if isempty(numericalIDs)
        disp('All participant data had already been converted from BIDS to .set');
        return;
    end

elseif ~strcmp(config.overwrite, 'on')
    warning('Undefined option for field config.overwrite - assume overwriting is on');
end

ftPath = fileparts(which('ft_defaults'));

% Import data set saved in BIDS, using a modified version of eeglab plugin
%--------------------------------------------------------------------------
try
    pop_editoptions('option_saveversion6', 0, 'option_single', 0, 'option_memmapdata', 0, 'option_savetwofiles', 1, 'option_storedisk', 0);
catch
    warning('Could NOT edit EEGLAB memory options!!');
end

pop_importbids_mobi(bidsDir,'datatypes',otherDataTypes,'outputdir', tempDir, 'participants', numericalIDs, 'matchchanlocs', config.match_electrodes_channels);

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
            sesMOTIONDir       = fullfile(tempDir, subjectDir, sesDir.name, 'motion');
            sesFilesEEG     = dir(sesEEGDir)';
            sesFilesBEH     = dir(sesBEHDir)';
            sesFilesMOTION  = dir(sesMOTIONDir)';
            allFiles        = [allFiles sesFilesEEG sesFilesBEH sesFilesMOTION];
        end

    else

        % for unisession, simply find all files in EEG folder
        eegDir          = fullfile(tempDir, subjectDir, 'eeg');
        behDir          = fullfile(tempDir, subjectDir, 'beh');
        motionDir       = fullfile(tempDir, subjectDir, 'motion');
        allFiles        = [dir(eegDir); dir(behDir); dir(motionDir)] ;

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

        if isMultiSession
            sessionName         = bidsNameSplit{strncmp(bidsNameSplit,'ses',3)}(5:end);
        end

        if find(strncmp(bidsNameSplit,'tracksys',8))
            isMultiTrackSys     = true;
            tracksysName        = bidsNameSplit{strncmp(bidsNameSplit,'tracksys',8)}(10:end);
        else
            isMultiTrackSys     = false;
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
                bemobilModality = 'MOTION';
                bidsFolderName = 'motion';
            case 'physio'
                bemobilModality = 'PHYSIO';
                bidsFolderName = 'beh';
            otherwise
                bemobilModality = bidsModality;
                disp(['Unknown modality' bidsModality ' saved as ' bemobilModality '.set'])
        end

        if isMultiSession

            dataDir         = fullfile(tempDir, subjectDir, ['ses-' sessionName] , bidsFolderName);

            % move files and then remove the empty eeg folder
            newDir          = fullfile(targetDir, [config.filename_prefix num2str(subjectNr)]);

            if ~isfolder(newDir)
                mkdir(newDir)
            end

            if strcmp(bemobilModality, 'MOTION') && isMultiTrackSys
                if isMultiRun
                    bemobilName     = [config.filename_prefix, num2str(subjectNr), '_' sessionName '_' bemobilModality, '_' tracksysName '_rec', runIndex, '_old', extension];
                else
                    bemobilName     = [config.filename_prefix, num2str(subjectNr), '_' sessionName '_' bemobilModality, '_' tracksysName '_old', extension];
                end
            else
                if isMultiRun
                    bemobilName     = [config.filename_prefix, num2str(subjectNr), '_' sessionName '_' bemobilModality, '_rec', runIndex, '_old', extension];
                else
                    bemobilName     = [config.filename_prefix, num2str(subjectNr), '_' sessionName '_' bemobilModality, '_old', extension];
                end
            end

            data    = pop_loadset('filepath', dataDir, 'filename', bidsName);
            pop_saveset(data, 'filepath', newDir, 'filename', bemobilName);

        else

            dataDir         = fullfile(tempDir, subjectDir, bidsFolderName);

            % move files and then remove the empty eeg folder
            newDir          = fullfile(targetDir, [config.filename_prefix num2str(subjectNr)]);

            if ~isfolder(newDir)
                mkdir(newDir)
            end

            % construct the new file name
            if strcmp(bemobilModality, 'MOTION') && isMultiTrackSys
                if isMultiRun
                    bemobilName     = [config.filename_prefix, num2str(subjectNr), '_' config.session_names{1} '_' bemobilModality '_' tracksysName  '_rec', runIndex '_old', extension];
                else
                    bemobilName     = [config.filename_prefix, num2str(subjectNr), '_' config.session_names{1} '_' bemobilModality '_' tracksysName '_old', extension];
                end
            else
                if isMultiRun
                    bemobilName     = [config.filename_prefix, num2str(subjectNr), '_' config.session_names{1} '_' bemobilModality, '_rec', runIndex '_old', extension];
                else
                    bemobilName     = [config.filename_prefix, num2str(subjectNr), '_' config.session_names{1} '_' bemobilModality '_old', extension];
                end
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
        subjectNr = str2double(nameArray{iN}(numel(config.filename_prefix) + 1:end));
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
    subjectNr = str2double(subDirList(iSub).name(numel(config.filename_prefix) + 1:end));

    % initialize a list of all tracking systems detected in data set
    trackingSystemsInData = {};

    % iterate over sessions
    for iSes = 1:numel(config.session_names)

        clear importfigs
        % find all EEG data
        eegFiles = {subjectFiles(contains({subjectFiles.name}, [config.session_names{iSes} '_EEG']) & contains({subjectFiles.name}, '_old.set')).name};
        eegFiles = natsortfiles(eegFiles);

        % resample and merge EEG
        %------------------------------------------------------------------
        if isempty(config.resample_freq)
            newSRate        = 250;
            warning('No resample frequency specified - data is still resampled to the default 250 Hz')
        else
            newSRate        = config.resample_freq;
        end

        % initializa a matrix storing EEG first and last time stamps (this is used to synch other streams)
        eegTimes = {};

        % initialize an empty cell array to store EEG events
        eegEvents = {};

        if numel(eegFiles) > 1

            % multi-run case
            EEGFileNameWithRun  = eegFiles{1};
            nameSplit           = regexp(EEGFileNameWithRun,'_', 'split'); % remove _rec entity
            nameJoined          = join(nameSplit(1:end-2),'_');
            EEGSessionFileName  = [nameJoined{1} '_old.set'];

            % loop over runs
            ALLEEG = []; CURRENTSET = [];
            for Ri = 1:numel(eegFiles)
                EEG                     = pop_loadset('filepath',fullfile(targetDir, subDirList(iSub).name),'filename', eegFiles{Ri});
                EEG.etc.effective_srate = EEG.srate;

                % set srate to nominal srate
                if isfield(EEG.etc,'nominal_srate')
                    srate_ratios(iSes,Ri) = EEG.etc.nominal_srate/EEG.srate;
                    disp(['Setting EEG srate from effective srate of ' num2str(EEG.srate) ' to nominal srate of ' num2str(EEG.etc.nominal_srate) ' Hz'])
                    assert(srate_ratios(iSes,1)>0.9 & srate_ratios(iSes,1)<1.1,'EEG effective and nominal srates are too different!')
                    EEG.srate = EEG.etc.nominal_srate;
                end
                
                % to prevent ft alt function from meddling with processing
                rmpath(genpath(fullfile(ftPath, [filesep 'external' filesep 'signal'])))
                
                % remove events with NaN latency before resampling
                eventLatencies      = [EEG.event(:).latency];
                nanInds             = find(isnan(eventLatencies));
                EEG.event(nanInds)  = []; 
                
                EEG                 = pop_resample( EEG, newSRate); % use filter-based resampling
                %                 [EEG]       = resampleToTime(EEG, newSRate, EEG.times(1), EEG.times(end), 0); % resample
                eegTimes{Ri}        = EEG.times;

                % round event times to have usable indices
                for i_event = 1:length(EEG.event)
                    EEG.event(i_event).latency = round(EEG.event(i_event).latency);
                end

                eegEvents{end +1}   = EEG.event;
                [ALLEEG,EEG,CURRENTSET]  = pop_newset(ALLEEG, EEG, CURRENTSET, 'study',0);

                % plot
                if ~isempty(EEG.event)
                    importfigs(Ri) = figure('color','w','position',[1 1 1920 1080]);

                    sgtitle(['Imported data from ' fullfile(EEG.filepath,erase(EEG.filename,'_old'))],'interpreter','none')

                    latencies_1 = EEG.event(1).latency;
                    latencies_2 = EEG.event(end).latency;

                    latenciesToPlot_1 = latencies_1-EEG.srate:latencies_1+2*EEG.srate;
                    latenciesToPlot_1(latenciesToPlot_1<1) = [];
                    latenciesToPlot_1(latenciesToPlot_1>EEG.pnts) = [];
                    latenciesToPlot_2 = latencies_2-EEG.srate:latencies_2+2*EEG.srate;
                    latenciesToPlot_2(latenciesToPlot_2<1) = [];
                    latenciesToPlot_2(latenciesToPlot_2>EEG.pnts) = [];

                    subplot(211); hold on; grid on; grid(gca,'minor')
                    title(['First event: "' EEG.event(1).type '"'],'interpreter','none')
                    yticks(-1)
                    yticklabels('')
                    xlabel('seconds')
                    xlim([EEG.times(latenciesToPlot_1(1)) EEG.times(latenciesToPlot_1(end))]/1000-EEG.times(latencies_1)/1000)

                    plot([0 0],[-1 100],'k')

                    my_yticks = yticks;
                    plot(EEG.times(latenciesToPlot_1)/1000-EEG.times(latencies_1)/1000,...
                        normalize(EEG.data(1,latenciesToPlot_1),...
                        'range',[my_yticks(end)+1 my_yticks(end)+2]), 'color', [78 165 216]/255)
                    yticks([yticks my_yticks(end)+1.5])
                    yticklabels([yticklabels
                        strrep(['EEG ' EEG.chanlocs(1).labels],'_', ' ')]);
                    ylim([-0.5 my_yticks(end)+2.5])

                    ax = gca;
                    ax.YAxis.MinorTickValues = ax.YAxis.Limits(1):0.2:ax.YAxis.Limits(2);

                    subplot(212); hold on; grid on; grid(gca,'minor')
                    title(['Last event: "' EEG.event(end).type '"'],'interpreter','none')
                    yticks(-1)
                    yticklabels('')
                    xlabel('seconds')
                    xlim([EEG.times(latenciesToPlot_2(1)) EEG.times(latenciesToPlot_2(end))]/1000-EEG.times(latencies_2)/1000)

                    plot([EEG.times(latencies_2) EEG.times(latencies_2)]/1000-EEG.times(latencies_2)/1000,[-1 100],'k')

                    my_yticks = yticks;
                    plot(EEG.times(latenciesToPlot_2)/1000-EEG.times(latencies_2)/1000,...
                        normalize(EEG.data(1,latenciesToPlot_2),...
                        'range',[my_yticks(end)+1 my_yticks(end)+2]), 'color', [78 165 216]/255)
                    yticks([yticks my_yticks(end)+1.5])
                    yticklabels([yticklabels
                        strrep(['EEG ' EEG.chanlocs(1).labels],'_', ' ')]);
                    ylim([-0.5 my_yticks(end)+2.5])

                    ax = gca;
                    ax.YAxis.MinorTickValues = ax.YAxis.Limits(1):0.2:ax.YAxis.Limits(2);
                else
                    warning('NO MARKERS IN THE FILE! NO PLOT CAN BE CREATED')
                end

            end
            [~, EEGMerged, ~]  = bemobil_merge(ALLEEG, EEG, CURRENTSET, 1:length(ALLEEG), EEGSessionFileName, fullfile(targetDir, [config.filename_prefix, num2str(subjectNr)]));
            EEG                = EEGMerged;

        elseif numel(eegFiles) == 1
            EEGSessionFileName      = eegFiles{1};
            EEG                     = pop_loadset('filepath',fullfile(targetDir, subDirList(iSub).name),'filename', eegFiles{1});
            EEG.etc.effective_srate = EEG.srate;

            % set srate to nominal srate
            if isfield(EEG.etc,'nominal_srate')
                srate_ratios(iSes,1) = EEG.etc.nominal_srate/EEG.srate;
                disp(['Setting EEG srate from effective srate of ' num2str(EEG.srate) ' to nominal srate of ' num2str(EEG.etc.nominal_srate) ' Hz'])
                assert(srate_ratios(iSes,1)>0.9 & srate_ratios(iSes,1)<1.1,'EEG effective and nominal srates are too different!')
                EEG.srate = EEG.etc.nominal_srate;
            end
            
            % to prevent ft alt function from meddling with processing
            rmpath(genpath(fullfile(ftPath, [filesep 'external' filesep 'signal'])))
                
            
            EEG                 = pop_resample( EEG, newSRate); % use filter-based resampling
            eegTimes            = EEG.times;
            
            % round event times to have usable indices
            for i_event = 1:length(EEG.event)
                EEG.event(i_event).latency = round(EEG.event(i_event).latency);
            end
            eegEvents{end +1}   = EEG.event;

            % plot
            if ~isempty(EEG.event)
                importfigs = figure('color','w','position',[1 1 1920 1080]);
                
                sgtitle(['Imported data from ' fullfile(EEG.filepath,erase(EEG.filename,'_old'))],'interpreter','none')
                
                latencies_1 = EEG.event(1).latency;
                latencies_2 = EEG.event(end).latency;
                
                latenciesToPlot_1 = latencies_1-EEG.srate:latencies_1+2*EEG.srate;
                latenciesToPlot_1(latenciesToPlot_1<1) = [];
                latenciesToPlot_1(latenciesToPlot_1>EEG.pnts) = [];
                latenciesToPlot_2 = latencies_2-EEG.srate:latencies_2+2*EEG.srate;
                latenciesToPlot_2(latenciesToPlot_2<1) = [];
                latenciesToPlot_2(latenciesToPlot_2>EEG.pnts) = [];
                
                subplot(211); hold on; grid on; grid(gca,'minor')
                title(['First event: "' EEG.event(1).type '"'],'interpreter','none')
                yticks(-1)
                yticklabels('')
                xlabel('seconds')
                xlim([EEG.times(latenciesToPlot_1(1)) EEG.times(latenciesToPlot_1(end))]/1000-EEG.times(latencies_1)/1000)
                
                plot([0 0],[-1 100],'k')
                
                my_yticks = yticks;
                plot(EEG.times(latenciesToPlot_1)/1000-EEG.times(latencies_1)/1000,...
                    normalize(EEG.data(1,latenciesToPlot_1),...
                    'range',[my_yticks(end)+1 my_yticks(end)+2]), 'color', [78 165 216]/255)
                yticks([yticks my_yticks(end)+1.5])
                yticklabels([yticklabels
                    strrep(['EEG ' EEG.chanlocs(1).labels],'_', ' ')]);
                ylim([-0.5 my_yticks(end)+2.5])
                
                ax = gca;
                ax.YAxis.MinorTickValues = ax.YAxis.Limits(1):0.2:ax.YAxis.Limits(2);
                
                subplot(212); hold on; grid on; grid(gca,'minor')
                title(['Last event: "' EEG.event(end).type '"'],'interpreter','none')
                yticks(-1)
                yticklabels('')
                xlabel('seconds')
                xlim([EEG.times(latenciesToPlot_2(1)) EEG.times(latenciesToPlot_2(end))]/1000-EEG.times(latencies_2)/1000)
                
                plot([EEG.times(latencies_2) EEG.times(latencies_2)]/1000-EEG.times(latencies_2)/1000,[-1 100],'k')
                
                my_yticks = yticks;
                plot(EEG.times(latenciesToPlot_2)/1000-EEG.times(latencies_2)/1000,...
                    normalize(EEG.data(1,latenciesToPlot_2),...
                    'range',[my_yticks(end)+1 my_yticks(end)+2]), 'color', [78 165 216]/255)
                yticks([yticks my_yticks(end)+1.5])
                yticklabels([yticklabels
                    strrep(['EEG ' EEG.chanlocs(1).labels],'_', ' ')]);
                ylim([-0.5 my_yticks(end)+2.5])
                
                ax = gca;
                ax.YAxis.MinorTickValues = ax.YAxis.Limits(1):0.2:ax.YAxis.Limits(2);
            else
                warning('NO MARKERS IN THE FILE! NO PLOT CAN BE CREATED')
            end

        else
            error(['No EEG file found in subject dir ' subDirList(iSub).name ', session ' config.session_names{iSes}] )
        end

        EEG = eeg_checkset(EEG, 'makeur');

        % save merged EEG file for the session
        EEG = pop_saveset(EEG, 'filename',[EEGSessionFileName(1:end-8) EEGSessionFileName(end-3:end)],'filepath',fullfile(targetDir, subDirList(iSub).name));
        disp(['Saved session file ' EEGSessionFileName(1:end-8) EEGSessionFileName(end-3:end)])


        % now iterate over other data types
        for iType = 1:numel(otherDataTypes)

            bemobilModality =  upper(otherDataTypes{iType});

            % resample and merge DATA
            %--------------------------------------------------------------
            if isempty(config.resample_freq)
                newSRate        = EEG.srate;
                warning(['No resample frequency specified -' bemobilModality 'data is still resampled to match EEG srate, ' num2str(newSRate) 'Hz'])
            else
                newSRate        = config.resample_freq;
            end

            % find all data of the type
            modalityFiles = {subjectFiles(contains({subjectFiles.name}, [config.session_names{iSes} '_' bemobilModality]) & contains({subjectFiles.name}, '_old.set')).name};

            trackingSystemsInSession = {};
            if strcmp(bemobilModality, 'MOTION') && isMultiTrackSys
                for MFi = 1:numel(modalityFiles)
                    % identify all tracking systems included in the session
                    modalityNamesSplit              = regexp(modalityFiles{MFi}, '_', 'split');
                    if contains(modalityNamesSplit{:,end-1}, 'rec')
                        trackingSystemsInSession{MFi}         = modalityNamesSplit{:,end-2};
                    else
                        trackingSystemsInSession{MFi}         = modalityNamesSplit{:,end-1};
                    end
                end
                trackingSystemsInSession = unique(trackingSystemsInSession);
                trackingSystemsInData   = unique([trackingSystemsInData trackingSystemsInSession]);
            else
                trackingSystemsInSession = {''};
            end

            for TSi = 1:numel(trackingSystemsInSession)

                dataFiles = modalityFiles(contains(modalityFiles, trackingSystemsInSession{TSi}));

                if numel(eegFiles) ~= numel(dataFiles)
                    warning(['Number of EEG files and data files of type ' bemobilModality ' do not match within session ''' config.session_names{iSes} ''''])
                end

                if numel(dataFiles) == 0
                    warning(['No data file of type ' bemobilModality ' found in session ''' config.session_names{iSes} '''. Skipping....'])
                    continue
                end

                if numel(dataFiles) > 1

                    DATAFileNameWithRun     = dataFiles{1};
                    nameSplit               = regexp(DATAFileNameWithRun,'_', 'split'); % remove _run entity
                    nameJoined              = join(nameSplit(1:end-2),'_');
                    DATASessionFileName      = [nameJoined{1} '_old.set'];

                    % loop over runs
                    ALLDATA = []; CURRENTSET = [];
                    for Ri = 1:numel(dataFiles)
                        DATA                        = pop_loadset('filepath',fullfile(targetDir, subDirList(iSub).name),'filename', dataFiles{Ri});
                        DATA.etc.effective_srate    = DATA.srate;
                        DATA                        = unwrapAngles(DATA); % unwrap angles before resampling
                        DATA.times                  = DATA.times/srate_ratios(iSes,Ri); % adjust for the slight imprecision in the EEG srate
                        if ~contains(lower(dataFiles{Ri}),lower(config.use_nominal_srate))

                            if DATA.srate > newSRate && config.do_aa_filter
                                disp('DATA srate is higher than new srate, performing anti-aliasing filter before downsampling.')
                                warning('IMPORTANT: The filter assumes equidistant samples, so if large sampling irregularities exist, this will be problematic! You can switch this function off with the "do_aa_filter" config entry (on by default).')
                                DATA = pop_eegfiltnew(DATA, [], newSRate/2*0.8); % low-pass filter such that the stopband edge is the neq nyquist freq
                            end

                            [DATA]      = resampleToTime(DATA, newSRate, eegTimes{Ri}, DATA.etc.starttime);
                        else
                            disp(['Data file ' dataFiles{Ri} config.use_nominal_srate ' using nominal srate!'])
                            assert(isfield(DATA.etc,'nominal_srate'),['Data file ' dataFiles{Ri} config.use_nominal_srate ' was specified to use nominal srate, but none was found!'])
                            DATA.srate  = DATA.etc.nominal_srate;
                            
                            % to prevent ft alt function from meddling with processing
                            rmpath(genpath(fullfile(ftPath, [filesep 'external' filesep 'signal'])))
                
                            DATA        = pop_resample( DATA, newSRate); % use filter-based resampling

                            % check beginning of data
                            sampleshift = round(DATA.etc.starttime*1000/(1/newSRate*1000));
                            if sampleshift < 0 % negative shift = need to cut data in the beginning
                                DATA.data   = DATA.data(:,1-sampleshift:end);
                            elseif sampleshift > 0 % positive shift = need to add nans
                                DATA.data   = [nan(DATA.nbchan,sampleshift) DATA.data];
                            end % shift of 0 means nothing needs to be changed

                            % check end of data
                            if size(DATA.data,2) > length(eegTimes{Ri}) % need to cut data short
                                DATA.data   = DATA.data(:,1:length(eegTimes{Ri}));
                            elseif size(DATA.data,2) < length(eegTimes{Ri}) % need to add nans
                                DATA.data   = [DATA.data nan(DATA.nbchan,length(eegTimes{Ri})-size(DATA.data,2))];
                            end % if length is perfect and shift is 0, nothing needs to be changed

                            % fix metadata
                            DATA.times  = eegTimes{Ri};
                            DATA.pnts   = length(eegTimes{Ri});
                            DATA.xmax   = DATA.times(end)/1000;
                            DATA.xmin   = 0;
                        end
                        DATA         = wrapAngles(DATA);
                        DATA.event   = eegEvents{Ri};
                        [ALLDATA,DATA,CURRENTSET]  = pop_newset(ALLDATA, DATA, CURRENTSET, 'study',0);
                        
                        % plot
                        if ~isempty(DATA.event)
                            figure(importfigs(Ri))
                            
                            latencies_1 = DATA.event(1).latency;
                            latencies_2 = DATA.event(end).latency;
                            
                            latenciesToPlot_1 = latencies_1-DATA.srate:latencies_1+2*DATA.srate;
                            latenciesToPlot_1(latenciesToPlot_1<1) = [];
                            latenciesToPlot_1(latenciesToPlot_1>DATA.pnts) = [];
                            latenciesToPlot_2 = latencies_2-DATA.srate:latencies_2+2*DATA.srate;
                            latenciesToPlot_2(latenciesToPlot_2<1) = [];
                            latenciesToPlot_2(latenciesToPlot_2>DATA.pnts) = [];
                            
                            idx = find(~contains({DATA.chanlocs.labels},'eul') & ~contains({DATA.chanlocs.labels},'quat') &...
                                ~contains({DATA.chanlocs.labels},'ori'),1,'first');
                            
                            subplot(211);
                            
                            my_yticks = yticks;
                            plot(DATA.times(latenciesToPlot_1)/1000-DATA.times(latencies_1)/1000,...
                                normalize(DATA.data(idx,latenciesToPlot_1),...
                                'range',[my_yticks(end)+1 my_yticks(end)+2]), 'color', [78 165 216]/255)
                            yticks([yticks my_yticks(end)+1.5])
                            yticklabels([yticklabels
                                strrep(['DATA ' DATA.chanlocs(idx).labels],'_', ' ')]);
                            ylim([-0.5 my_yticks(end)+2.5])
                            
                            ax = gca;
                            ax.YAxis.MinorTickValues = ax.YAxis.Limits(1):0.2:ax.YAxis.Limits(2);
                            
                            subplot(212);
                            
                            my_yticks = yticks;
                            plot(DATA.times(latenciesToPlot_2)/1000-DATA.times(latencies_2)/1000,...
                                normalize(DATA.data(idx,latenciesToPlot_2),...
                                'range',[my_yticks(end)+1 my_yticks(end)+2]), 'color', [78 165 216]/255)
                            yticks([yticks my_yticks(end)+1.5])
                            yticklabels([yticklabels
                                strrep(['DATA ' DATA.chanlocs(idx).labels],'_', ' ')]);
                            ylim([-0.5 my_yticks(end)+2.5])
                            
                            ax = gca;
                            ax.YAxis.MinorTickValues = ax.YAxis.Limits(1):0.2:ax.YAxis.Limits(2);
                        else
                            warning('NO MARKERS IN THE FILE! NO PLOT CAN BE CREATED')
                        end
                            

                    end
                    [~, DATAMerged, ~]  = bemobil_merge(ALLDATA, DATA, CURRENTSET, 1:length(ALLDATA), DATASessionFileName, fullfile(targetDir, [config.filename_prefix, num2str(subjectNr)]));
                    DATA             = DATAMerged;
                elseif numel(dataFiles) == 1
                    DATASessionFileName  = dataFiles{1};
                    DATA                        = pop_loadset('filepath', fullfile(targetDir, subDirList(iSub).name), 'filename', dataFiles{1});
                    DATA.etc.effective_srate    = DATA.srate;
                    DATA                        = unwrapAngles(DATA);
                    DATA.times                  = DATA.times/srate_ratios(iSes); % adjust for the slight imprecision in the EEG srate
                    if ~contains(lower(dataFiles{1}),lower(config.use_nominal_srate))

                        if DATA.srate > newSRate && config.do_aa_filter
                            disp('DATA srate is higher than new srate, performing anti-aliasing filter before downsampling.')
                            warning('IMPORTANT: The filter assumes equidistant samples, so if large sampling irregularities exist, this will be problematic! You can switch this function off with the "do_aa_filter" config entry (on by default).')
                            DATA = pop_eegfiltnew(DATA, [], newSRate/2*0.8); % low-pass filter such that the stopband edge is the neq nyquist freq
                        end

                        [DATA]       = resampleToTime(DATA, newSRate, eegTimes, DATA.etc.starttime);
                    else
                        disp(['Data file ' dataFiles{1} ' using nominal srate!'])
                        assert(isfield(DATA.etc,'nominal_srate'),['Data file ' dataFiles{1} ' was specified to use nominal srate, but none was found!'])
                        DATA.srate  = DATA.etc.nominal_srate;
                        
                        % to prevent ft alt function from meddling with processing
                        rmpath(genpath(fullfile(ftPath, [filesep 'external' filesep 'signal'])))

                        DATA        = pop_resample( DATA, newSRate); % use filter-based resampling

                        % check beginning of data
                        sampleshift = round(DATA.etc.starttime*1000/(1/newSRate*1000));
                        if sampleshift < 0 % negative shift = need to cut data in the beginning
                            DATA.data   = DATA.data(:,1-sampleshift:end);
                        elseif sampleshift > 0 % positive shift = need to add nans
                            DATA.data   = [nan(DATA.nbchan,sampleshift) DATA.data];
                        end % shift of 0 means nothing needs to be changed

                        % check end of data
                        if size(DATA.data,2) > length(eegTimes) % need to cut data short
                            DATA.data   = DATA.data(:,1:length(eegTimes));
                        elseif size(DATA.data,2) < length(eegTimes) % need to add nans
                            DATA.data   = [DATA.data nan(DATA.nbchan,length(eegTimes)-size(DATA.data,2))];
                        end % if length is perfect and shift is 0, nothing needs to be changed

                        % fix metadata
                        DATA.times  = eegTimes;
                        DATA.pnts   = length(eegTimes);
                        DATA.xmax   = DATA.times(end)/1000;
                        DATA.xmin   = 0;
                    end
                    DATA             = wrapAngles(DATA);
                    DATA.event       = eegEvents{1};
                    
                    % plot
                    if ~isempty(DATA.event)
                        figure(importfigs)
                        
                        latencies_1 = DATA.event(1).latency;
                        latencies_2 = DATA.event(end).latency;
                        
                        latenciesToPlot_1 = latencies_1-DATA.srate:latencies_1+2*DATA.srate;
                        latenciesToPlot_1(latenciesToPlot_1<1) = [];
                        latenciesToPlot_1(latenciesToPlot_1>DATA.pnts) = [];
                        latenciesToPlot_2 = latencies_2-DATA.srate:latencies_2+2*DATA.srate;
                        latenciesToPlot_2(latenciesToPlot_2<1) = [];
                        latenciesToPlot_2(latenciesToPlot_2>DATA.pnts) = [];
                        
                        idx = find(~contains({DATA.chanlocs.labels},'eul') & ~contains({DATA.chanlocs.labels},'quat') &...
                            ~contains({DATA.chanlocs.labels},'ori'),1,'first');
                        
                        subplot(211);
                        
                        my_yticks = yticks;
                        plot(DATA.times(latenciesToPlot_1)/1000-DATA.times(latencies_1)/1000,...
                            normalize(DATA.data(idx,latenciesToPlot_1),...
                            'range',[my_yticks(end)+1 my_yticks(end)+2]), 'color', [78 165 216]/255)
                        yticks([yticks my_yticks(end)+1.5])
                        yticklabels([yticklabels
                            strrep(['DATA ' DATA.chanlocs(idx).labels],'_', ' ')]);
                        ylim([-0.5 my_yticks(end)+2.5])
                        
                        ax = gca;
                        ax.YAxis.MinorTickValues = ax.YAxis.Limits(1):0.2:ax.YAxis.Limits(2);
                        
                        subplot(212);
                        
                        my_yticks = yticks;
                        plot(DATA.times(latenciesToPlot_2)/1000-DATA.times(latencies_2)/1000,...
                            normalize(DATA.data(idx,latenciesToPlot_2),...
                            'range',[my_yticks(end)+1 my_yticks(end)+2]), 'color', [78 165 216]/255)
                        yticks([yticks my_yticks(end)+1.5])
                        yticklabels([yticklabels
                            strrep(['DATA ' DATA.chanlocs(idx).labels],'_', ' ')]);
                        ylim([-0.5 my_yticks(end)+2.5])
                        
                        ax = gca;
                        ax.YAxis.MinorTickValues = ax.YAxis.Limits(1):0.2:ax.YAxis.Limits(2);
                    else
                        warning('NO MARKERS IN THE FILE! NO PLOT CAN BE CREATED')
                    end

                else
                    warning(['No file of modality ' bemobilModality ' found in subject dir ' subDirList(iSub).name ', session ' config.session_names{iSes}] )
                end

                % create urevents
                DATA             = eeg_checkset(DATA, 'makeur');

                % save merged data file for the session
                DATA = pop_saveset(DATA, 'filename', [DATASessionFileName(1:end-8) DATASessionFileName(end-3:end)],'filepath',fullfile(targetDir, subDirList(iSub).name));
                disp(['Saved session file ' [DATASessionFileName(1:end-8) DATASessionFileName(end-3:end)]])

            end
        end

        % save fig
        if exist('importfigs','var')
            if length(importfigs) == 1
                savefig(importfigs,...
                    fullfile(targetDir, subDirList(iSub).name, [subDirList(iSub).name '_' config.session_names{iSes} '_imported-data']))
                print(importfigs,...
                    fullfile(targetDir, subDirList(iSub).name, [subDirList(iSub).name '_' config.session_names{iSes} '_imported-data']),'-dpng')
                close(importfigs)
            else
                for i = 1:length(importfigs)
                    if ishandle(importfigs(i))
                        savefig(importfigs(i),...
                            fullfile(targetDir, subDirList(iSub).name, [subDirList(iSub).name '_' config.session_names{iSes} '_run-' num2str(i) '_imported-data']))
                        print(importfigs(i),...
                            fullfile(targetDir, subDirList(iSub).name, [subDirList(iSub).name '_' config.session_names{iSes} '_run-' num2str(i) '_imported-data']),'-dpng')
                        close(importfigs(i))
                    end
                end
            end
        end

    end

    % remove unnecessary files prior to merging
    subjectFiles = dir(fullfile(targetDir, subDirList(iSub).name));
    toDelete = {subjectFiles(contains({subjectFiles.name}, '_old.fdt') | contains({subjectFiles.name}, '_old.set')).name};

    for iD = 1:numel(toDelete)
        delete(fullfile(targetDir, subDirList(iSub).name, toDelete{iD}));
    end

    % merge sessions
    %----------------------------------------------------------------------
    if isMultiSession
        if numel(config.session_names) > 1

            % EEG data
            %--------------------------------------------------------------
            ALLEEG = []; CURRENTSET = [];
            for Si = 1:numel(config.session_names)
                [outPath, outName] = sessionfilename(targetDir,'EEG', config, Si, subjectNr);
                EEG         = pop_loadset('filepath',outPath,'filename',outName);
                [ALLEEG,EEG,CURRENTSET]  = pop_newset(ALLEEG, EEG, CURRENTSET, 'study',0);
            end
            [~, ~, ~]  = bemobil_merge(ALLEEG, EEG, CURRENTSET, 1:length(ALLEEG), [config.filename_prefix, num2str(subjectNr), '_' config.merged_filename '_EEG.set'], fullfile(targetDir, [config.filename_prefix, num2str(subjectNr)]));

            % MOTION and other data
            %--------------------------------------------------------------
            motionsets = {};
            motionsets_tracksys = {};
            for iType = 1:numel(otherDataTypes)

                if isempty(trackingSystemsInData)
                    trackingSystemsInData = {''};
                end

                for TSi = 1:numel(trackingSystemsInData)

                    fillIndices = zeros(1,numel(config.session_names));

                    for Si = 1:numel(config.session_names)
                        if strcmpi(otherDataTypes{iType}, 'MOTION')
                            [outPath, outName] = sessionfilename(targetDir,['MOTION_' trackingSystemsInData{TSi}], config, Si, subjectNr);
                        else
                            [outPath, outName] = sessionfilename(targetDir, upper(otherDataTypes{iType}), config, Si, subjectNr);
                        end

                        % check existence of the session file
                        if exist(fullfile(outPath, outName), 'file')
                            DATA         = pop_loadset('filepath',outPath,'filename',outName);
                        else
                            fillIndices(Si) = 1;
                            disp(['File ' outName ' does not exist - filling with NaNs for merging over sessions']);
                        end
                    end

                    % fill with NaNs for merging across sessions, if necessary
                    fillDataArray = {};
                    if isempty(DATA)
                        if isempty(trackingSystemsInData{1})
                            warning(['No ' upper(otherDataTypes{iType}) ' files found for any session, unable to merge over sessions'])
                        else
                            warning(['No ' upper(otherDataTypes{iType}) ' files in ' trackingSystemsInData{TSi} 'system found for any session, unable to merge'])
                        end
                    else
                        for Si = 1:numel(config.session_names)
                            if fillIndices(Si) % if data needs filling
                                % load again the EEG data in order to check the length and srate
                                [outPath, outName] = sessionfilename(targetDir,'EEG', config, Si, subjectNr);
                                fillData    = pop_loadset('filepath',outPath,'filename',outName);

                                % take any loaded data in order to check the channel information
                                fillData.nbchan     = DATA.nbchan;
                                fillData.data       = NaN(fillData.nbchan, size(fillData.data,2));
                                fillData.chanlocs   = DATA.chanlocs;
                                fillData.setname    = [];
                                fillData.filename   = [];
                                fillData            = eeg_checkset(fillData);
                                fillDataArray{Si}   = fillData;
                            end
                        end
                    end

                    ALLDATA = []; CURRENTSET = []; DATA = [];
                    for Si = 1:numel(config.session_names)
                        if strcmpi(otherDataTypes{iType}, 'MOTION')
                            [outPath, outName] = sessionfilename(targetDir,['MOTION_' trackingSystemsInData{TSi}], config, Si, subjectNr);
                        else
                            [outPath, outName] = sessionfilename(targetDir, upper(otherDataTypes{iType}), config, Si, subjectNr);
                        end

                        % check existence of the session file
                        if exist(fullfile(outPath, outName), 'file')
                            DATA         = pop_loadset('filepath',outPath,'filename',outName);
                            [ALLDATA,DATA,CURRENTSET]  = pop_newset(ALLDATA, DATA, CURRENTSET, 'study',0);
                        else
                            DATA        = fillDataArray{Si};
                            [ALLDATA,DATA,CURRENTSET]  = pop_newset(ALLDATA, DATA, CURRENTSET, 'study',0);
                        end
                    end
                    if strcmpi(otherDataTypes{iType}, 'MOTION')
                        % store motion sets separately and merge all into one below
                        [~, motionsets{end+1}, ~]  = bemobil_merge(ALLDATA, DATA, CURRENTSET, 1:length(ALLDATA));
                        motionsets_tracksys{end+1} = trackingSystemsInData{TSi};
                    else
                        [~, ~, ~]  = bemobil_merge(ALLDATA, DATA, CURRENTSET, 1:length(ALLDATA), [config.filename_prefix, num2str(subjectNr), '_' config.merged_filename '_' upper(otherDataTypes{iType}) '.set'], fullfile(targetDir, [config.filename_prefix, num2str(subjectNr)]));
                    end
                end
            end
            if ~isempty(motionsets)
                % motion sets were found and now need to be merged over tracking systems

                [~, sorting_idx] = sort(motionsets_tracksys);
                motionsets = motionsets(sorting_idx);

                mergedmotion = motionsets{1};

                for i = 2:length(motionsets)
                    mergedmotion.nbchan = mergedmotion.nbchan + motionsets{i}.nbchan;
                    mergedmotion.chanlocs = [mergedmotion.chanlocs motionsets{i}.chanlocs];
                    mergedmotion.data = [mergedmotion.data; motionsets{i}.data];
                end

                pop_saveset(mergedmotion, 'filepath', fullfile(targetDir, [config.filename_prefix, num2str(subjectNr)]),...
                    'filename', [config.filename_prefix, num2str(subjectNr), '_' config.merged_filename '_MOTION.set']);

            end
        else
            % still merge all MOTION datasets
            if ~isempty(trackingSystemsInData)

                motionsets = {};
                motionsets_tracksys = {};

                for TSi = 1:numel(trackingSystemsInData)
                    [outPath, outName] = sessionfilename(targetDir,['MOTION_' trackingSystemsInData{TSi}], config, 1, subjectNr);

                    % load the tracksys file
                    if exist(fullfile(outPath, outName), 'file')
                        motionsets{end+1}         = pop_loadset('filepath',outPath,'filename',outName);
                    else
                        disp(['File ' outName ' does not exist - this error can not happen, what did you do? Please check in with the BeMoBIL...']);
                    end

                    motionsets_tracksys{end+1} = trackingSystemsInData{TSi};
                end
                if ~isempty(motionsets)
                    % motion sets were found and now need to be merged over tracking systems

                    [~, sorting_idx] = sort(motionsets_tracksys);
                    motionsets = motionsets(sorting_idx);

                    mergedmotion = motionsets{1};

                    for i = 2:length(motionsets)
                        mergedmotion.nbchan = mergedmotion.nbchan + motionsets{i}.nbchan;
                        mergedmotion.chanlocs = [mergedmotion.chanlocs motionsets{i}.chanlocs];
                        mergedmotion.data = [mergedmotion.data; motionsets{i}.data];
                    end

                    [outPath, outName] = sessionfilename(targetDir,['MOTION'], config, 1, subjectNr);

                    pop_saveset(mergedmotion, 'filepath', outPath, 'filename', outName);

                end
            end
        end
    else
        % still merge all MOTION datasets
        if ~isempty(trackingSystemsInData)

            motionsets = {};
            motionsets_tracksys = {};

            for TSi = 1:numel(trackingSystemsInData)
                [outPath, outName] = sessionfilename(targetDir,['MOTION_' trackingSystemsInData{TSi}], config, 1, subjectNr);

                % load the tracksys file
                if exist(fullfile(outPath, outName), 'file')
                    motionsets{end+1}         = pop_loadset('filepath',outPath,'filename',outName);
                else
                    disp(['File ' outName ' does not exist - this error can not happen, what did you do? Please check in with the BeMoBIL...']);
                end

                motionsets_tracksys{end+1} = trackingSystemsInData{TSi};
            end
            if ~isempty(motionsets)
                % motion sets were found and now need to be merged over tracking systems

                [~, sorting_idx] = sort(motionsets_tracksys);
                motionsets = motionsets(sorting_idx);

                mergedmotion = motionsets{1};

                for i = 2:length(motionsets)
                    mergedmotion.nbchan = mergedmotion.nbchan + motionsets{i}.nbchan;
                    mergedmotion.chanlocs = [mergedmotion.chanlocs motionsets{i}.chanlocs];
                    mergedmotion.data = [mergedmotion.data; motionsets{i}.data];
                end

                [outPath, outName] = sessionfilename(targetDir,['MOTION'], config, 1, subjectNr);

                pop_saveset(mergedmotion, 'filepath', outPath, 'filename', outName);

            end
        end
    end
end

disp('BIDS to .set conversion finished')

end

function [outPath, outName] = sessionfilename(targetDir, modality, bemobil_config, sesnr, subnr)

outName     = [bemobil_config.filename_prefix, num2str(subnr), '_', bemobil_config.session_names{sesnr} '_' modality '.set'];
outPath     = fullfile(targetDir,[bemobil_config.filename_prefix, num2str(subnr)]);

end

function [outEEG] = resampleToTime(EEG, newSRate, newTimes, offset)
% offset is in seconds
%--------------------------------------------------------------------------

% save old times
oldTimes                = EEG.times;

% save old nan beginning and end
nanbegin = oldTimes(find(~any(isnan(EEG.data)),1,'first'));
nanend = oldTimes(find(~any(isnan(EEG.data)),1,'last'));

% Note that in fieldtrip time is in seconds
newTimesFT                = newTimes/1000;

resamplecfg.time        = {newTimesFT};
resamplecfg.detrend     = 'no';
resamplecfg.method      = 'pchip';
resamplecfg.extrapval   = nan;
EEG.group = 1; EEG.condition = 1;
ftData                  = eeglab2fieldtrip( EEG, 'raw', 'none' );
ftData.time{1}          = oldTimes/1000 + offset;
resampledData           = ft_resampledata(resamplecfg, ftData);
EEG.data                = resampledData.trial{1};
EEG.srate               = newSRate;
EEG.times               = newTimes; % convert back to miliseconds
EEG.pnts                = size(EEG.data,2);
EEG.urevent             = EEG.event;

% sometimes broken recordings lead to NaN values in the latency fields
for iE = 1:numel(EEG.event)
    if isnan(EEG.urevent(iE).latency)
        warning(['NaN value detected in event ' num2str(iE)])
    else
        EEG.event(iE).latency        = find(EEG.times > oldTimes(round(EEG.urevent(iE).latency)),1,'first');
    end
end

EEG.setname = EEG.filename(1:end-8);

% replace extrapolated values with nan (it shouldnt do it but for some reason it does...)
EEG.data(:,1:find(EEG.times>nanbegin+offset*1000,1,'first')-1) = deal(nan); % offset added because its also added to the interpolation
EEG.data(:,find(EEG.times<nanend+offset*1000,1,'last')+1:end) = deal(nan);

% checkset
outEEG = eeg_checkset(EEG, 'makeur');

end

function [DATA] = unwrapAngles(DATA)

% unwrap any kind of angular data before resampling
angleind = [];
for Ci = 1:numel(DATA.chanlocs)
    if ischar(DATA.chanlocs(Ci).units)
        if  contains(lower(DATA.chanlocs(Ci).units), 'rad') || contains(lower(DATA.chanlocs(Ci).units), 'deg')
            angleind = [angleind Ci];
        end
    end
end

DATA.data(angleind,:)   = unwrap(DATA.data(angleind,:),[],2);

end


function [DATA] = wrapAngles(DATA)

% wrap any kind of angular data before resampling
angleind = [];
for Ci = 1:numel(DATA.chanlocs)
    if ~isempty(DATA.chanlocs(Ci).units)
        if  contains(lower(DATA.chanlocs(Ci).units), 'rad') || contains(lower(DATA.chanlocs(Ci).units), 'deg')
            angleind = [angleind Ci];
        end
    end
end

DATA.data(angleind,:)   = wrapToPi(DATA.data(angleind,:));

end


%--------------------------------------------------------------------------
function [newconfig] =  checkfield(oldconfig, fieldName, defaultValue, defaultValueText)
% This function checks for existence of required fields
% for optional fields that have default values, assign them

% throw an error if a required field has not been specified
if strcmp(defaultValue, 'required')
    if ~isfield(oldconfig, fieldName)
        error(['Required config field ' fieldName ' not specified'])
    end
end

% assign default value to an optional field
newconfig   = oldconfig;

if ~isfield(oldconfig, fieldName)
    newconfig.(fieldName) = defaultValue;
    warning(['Config field ' fieldName ' not specified- using default value: ' defaultValueText])
end


end
