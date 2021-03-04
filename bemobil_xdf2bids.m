function bemobil_xdf2bids(bemobil_config, numericalIDs)

% Examaple script for converting .xdf recordings to BIDS
% see bemobil_bidstools.md 
% 
% Inputs : 
%   bemobil_config
%   example fields 
%       bemobil_config.study_folder             = 'E:\Project_BIDS\example_dataset_MWM\';
%       bemobil_config.filename_prefix          = 'sub_';
%       bemobil_config.raw_data_folder          = '0_raw-data\';
%       bemobil_config.bids_data_folder         = '1_BIDS-data\'; 
%       bemobil_config.raw_EEGLAB_data_folder   = '2_basic-EEGLAB\';
%       bemobil_config.filenames                = {'VR' 'desktop'}; 
%       bemobil_config.rigidbody_streams        = {'playerTransform','playerTransfom','rightHand', 'leftHand', 'Torso'};
%       bemobil_config.bids_rbsessions          = [1,1; 1,1; 0,1; 0,1; 0,1]; 
%                                                  logicals : indicate which rbstreams are present in which sessions
%       bemobil_config.bids_eegkeyword          = {'BrainVision'};
%                                                  a unique keyword used to identify the eeg stream in the .xdf file             
%       bemobil_config.bids_tasklabel           = 'VNE1';
%       bemobil_config.channel_locations_filename = 'VN_E1_eloc.elc'; 
%   
%   numericalIDs
%       array of participant numerical IDs in the data set
%  
%   Note on multi-session and multi-run data sets :  
%
%           Entries in bemobil_config.filenames will search through the 
%           raw data directory of the participant and group together
%           .xdf files with matching keyword in the name into one session
%           If there are multiple files in one session, they will be given
%           separate 'run' numbers in the file name
%           !! IMPORTANT !!
%           The order of runs rely on incremental name sorting 
%           using the "Nature Order File Sorting Tool" - Stephen Cobeldick (2021). Natural-Order Filename Sort (https://www.mathworks.com/matlabcentral/fileexchange/47434-natural-order-filename-sort), MATLAB Central File Exchange. Retrieved March 1, 2021.
%
%           for example,
%                   sub-1\sub-1_VNE1_VR_rec1.xdf
%                   sub-1\sub-1_VNE1_VR_rec2.xdf
%                   sub-1\sub-1_VNE1_desktop.xdf
%           will be organized into
%                   sub-001\ses-VR\sub-001_ses-VR_task-VNE1_run-1_eeg.bdf
%                   sub-001\ses-VR\sub-001_ses-VR_task-VNE1_run-2_eeg.bdf
%                   sub-001\ses-desktop\sub-001_ses-desktop_task-VNE1_eeg.bdf
%
% Author : Sein Jeung (seinjeung@gmail.com)
%--------------------------------------------------------------------------

% input check and default value assignment 
%--------------------------------------------------------------------------
if ~isfield(bemobil_config, 'bids_data_folder')
    bemobil_config.bids_data_folder = '1_BIDS-data\';
end

if ~isfield(bemobil_config, 'bids_rbsessions')
    bemobil_config.bids_rbsessions    = true(numel(bemobil_config.filenames),numel(bemobil_config.rigidbody_streams)); 
else
    if ~islogical(bemobil_config.bids_rbsessions)
        bemobil_config.bids_rbsessions = logical(bemobil_config.bids_rbsessions); 
    end
end

if ~isfield(bemobil_config, 'bids_eegkeyword')
    bemobil_config.bids_rbsessions      = {'BrainVision'}; 
end

if ~isfield(bemobil_config, 'bids_tasklabel')
    bemobil_config.bids_tasklabel       = 'defaulttask';
end


% start processing 
%--------------------------------------------------------------------------
% initialize fieldtrip
ft_defaults

% add natsortfiles to path
[filepath,~,~] = fileparts(which('bemobil_xdf2bids')); 
addpath(fullfile(filepath, 'resources', 'natsortfiles'))


% path to sourcedata
sourceDataPath                          = fullfile(bemobil_config.study_folder, bemobil_config.raw_data_folder(1:end-1));
addpath(genpath(sourceDataPath))

% range of (effective) sampling rate of data, used for identifying streams
% in .xdf
motionStreamNames                       = bemobil_config.rigidbody_streams;
eegStreamName                           = bemobil_config.bids_eegkeyword;

% funcions that resolve dataset-specific problems
motionCustom                            = 'bemobil_bidsmotion';

% load general metadata that apply to all participants
bemobil_bidsconfig_general;

%--------------------------------------------------------------------------
% loop over participants
for pi = 1:numel(numericalIDs)
    
    
    participantNr   = numericalIDs(pi); 
    participantDir  = fullfile(sourceDataPath, [bemobil_config.filename_prefix num2str(participantNr)]);
    
    % find all .xdf files for the given session in the participant directory
    participantFiles    = dir(participantDir);
    fileNameArray       = {participantFiles.name};
    
    % loop over sessions
    for si = 1:numel(bemobil_config.filenames)
        
        sessionFiles        = participantFiles(contains(fileNameArray, '.xdf') & contains(fileNameArray, bemobil_config.filenames{si}));
        
        % sort files by natural order 
        sortedFileNames     = natsortfiles({sessionFiles.name}); 
        
        % loop over files in each session. 
        % Here 'di' will index files as runs.
        for di = 1:numel(sortedFileNames)
            
            % construct file and participant- and file- specific config
            bemobil_bidsconfig_participants;
            
            %--------------------------------------------------------------
            %                Convert Motion Data to BIDS
            %--------------------------------------------------------------
            % import motion data
            motionSource                = xdf2fieldtrip(cfg.dataset,'streamkeywords', motionStreamNames(bemobil_config.bids_rbsessions(si,:)));
            
            % if needed, execute a custom function for any alteration
            % to the data to address dataset specific issues
            % (quat2eul conversion, for instance)
            if exist('motionCustom','var')
                motion                      = feval(motionCustom, motionSource, motionStreamNames(bemobil_config.bids_rbsessions(si,:)));
            else
                motion = motionSource;
            end
            
            % construct motion metadata
            bemobil_bidsconfig_motion;
            
            % write motion files in bids format
            data2bids(motioncfg, motion);
            
            %--------------------------------------------------------------
            %                  Convert EEG Data to BIDS
            %--------------------------------------------------------------
            % import eeg data
            eegSource                       = xdf2fieldtrip(cfg.dataset,'streamkeywords', eegStreamName);
            
            % if needed, execute a custom function for any alteration
            % to the data to address dataset specific issues
            if exist('eegCustom','var')
                eeg                      = feval(eegCustom, eegSource);
            else
                eeg = eegSource;
            end
            
            % construct eeg metadata
            bemobil_bidsconfig_eeg;
            
            % read in the event stream (synched to the EEG stream)
            eegcfg.event = ft_read_event(cfg.dataset);
            
            % write eeg files in bids format
            data2bids(eegcfg, eeg);
        end
    end
end


end