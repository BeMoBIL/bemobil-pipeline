% Process a dataset containing rigid body data streams (which contain 6 DOF data of position and orientation).
% Orientation can be either euler angles or unit quaternions. Orientation channels need to have the suffixes 'eul_x/y/z'
% or 'quat_x/y/z/w'. The dataset must contain information about the rigid bodies in the chanlocs substruct, specifically
% the respective rigid bodies need to be marked by the same "EEG.chanlocs.tracked_location" entry. Processing contains
% optional lowpass filtering, transformation to euler angles, and taking the first two derivatives (velocity and
% acceleration).
%
% Inputs:
%   ALLEEG                    - complete EEGLAB data set structure
%   EEG_motion_in             - EEGLAB EEG structure that should be processed. Best to have all blocks merged into one
%                                file.
%   CURRENTSET                - index of current EEGLAB EEG structure within ALLEEG
%   subject                   - subject number of the current subject (necessary for filepaths and storage)
%   bemobil_config            - configuration struct with all necessary information. See EEG_processing_example file
%                                that comes with this function!
%
% Outputs:
%   ALLEEG                    - complete EEGLAB data set structure
%   EEG_motion_out            - current EEGLAB EEG structure
%   Currentset                - index of current EEGLAB EEG structure within ALLEEG
%
%   .set data file of current EEGLAB EEG structure stored on disk
%
% See also:
%   EEGLAB
%
% Authors: Marius Klug, 2021

function [ALLEEG, EEG_motion_out, CURRENTSET] = bemobil_process_all_motion(ALLEEG, EEG_motion_in, CURRENTSET,...
    subject, bemobil_config, force_recompute)

% check config
bemobil_config = bemobil_check_config(bemobil_config);

% make sure the data is stored in double precision, large datafiles are supported, no memory mapped objects are
% used but data is processed locally, and two files are used for storing sets (.set and .fdt)
try 
    pop_editoptions('option_saveversion6', 0, 'option_single', 0, 'option_memmapdata', 0, 'option_savetwofiles', 1, 'option_storedisk', 0);
catch
    warning('Could NOT edit EEGLAB memory options!!'); 
end

if ~exist('force_recompute','var')
    force_recompute = false;
end
if force_recompute
    warning('RECOMPUTING OLD FILES IF FOUND!!!')
end


%% check if the entire processing was done already
output_filepath = fullfile(bemobil_config.study_folder, bemobil_config.motion_analysis_folder, [bemobil_config.filename_prefix num2str(subject)]);

if ~force_recompute
    try
        
        EEG = pop_loadset('filename', [bemobil_config.filename_prefix num2str(subject) '_'...
            bemobil_config.processed_motion_filename], 'filepath', output_filepath);
        [ALLEEG, EEG_motion_out, CURRENTSET] = pop_newset(ALLEEG, EEG, 0,'study',0);
        
        warning('Old processed motion file already existed, skipping this processing!')
        
    catch
        disp('...failed. Computing now.')
    end
end


if ~exist('EEG_motion_out','var')
    
    %%
    all_tracked_points = lower({EEG_motion_in.chanlocs.tracked_point}');
    all_rigidbody = unique(all_tracked_points);
    
    EEG_motion_out = EEG_motion_in;
    EEG_motion_out.nbchan = 0;
    EEG_motion_out.chanlocs = [];
    EEG_motion_out.data = [];
    
    for i_rigidbody = 1:length(all_rigidbody)
        
        disp('-----------------------------------------------------------------')
        disp([num2str(i_rigidbody) '/' num2str(length(all_rigidbody)) ': Processing ' all_rigidbody{i_rigidbody} ' rigid body data!'])
        disp('-----------------------------------------------------------------')
        
        idx_this_rb = find(contains(all_tracked_points,all_rigidbody{i_rigidbody}));
%         assert(length(idx_this_rb)==6||length(idx_this_rb)==7,'Rigidbody does not contain 6 (Euler) or 7 (quaternions) channels!')
        EEG_single_rbs = pop_select( EEG_motion_in, 'channel',idx_this_rb);
        
        EEG_single_rbs = bemobil_motion_process_rigidbody_data(EEG_single_rbs,bemobil_config.lowpass_motion,bemobil_config.lowpass_motion_after_derivative);
        
        EEG_motion_out.nbchan = EEG_motion_out.nbchan + EEG_single_rbs.nbchan;
        EEG_motion_out.chanlocs = [EEG_motion_out.chanlocs EEG_single_rbs.chanlocs];
        EEG_motion_out.data = [EEG_motion_out.data; EEG_single_rbs.data];
    end
    
    %%
    mkdir(output_filepath)
    pop_saveset(EEG_motion_out, 'filename', [ bemobil_config.filename_prefix num2str(subject)...
        '_' bemobil_config.processed_motion_filename], 'filepath', output_filepath);
    
end
