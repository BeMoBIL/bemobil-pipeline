function EEG_mocap_out = bemobil_process_mocap_data(EEG_mocap_in,lowpass_mocap,lowpass_after_derivative)

% Process a single rigid body data stream which contains 6 DOF data of position and orientation. Orientation can be
% either euler angles. Channels need to have the suffixes 'euler_x/y/z'), or unit quaternions (channels need to have the
% suffixes 'quat_x/y/z/w'. Processing contains lowpass filtering, transformation to euler angles, and taking the first
% two derivatives (velocity and acceleration).
%
% Input arguments:
%       EEG_mocap_in                - EEGLAB dataset containing mocap channels with 6 DOF
%       lowpass_mocap               - lowpass passband edge for filtering the raw mocap data. If empty, no filter is 
%                                       applied. 7Hz is recommended.
%       lowpass_after_derivative    - lowpass passband edge for filtering during after time derivative. If empty, no 
%                                       filter is applied. Can help if the derivatives are too noisy, but not
%                                       specifically recommended.
%
% Output argument:
%       EEG_mocap_out   - EEGLAB dataset containing processed mocap channels with 6 DOF in euler angles, plus the first
%       two derivatives (18 channels altogether)

all_tracked_points = lower({EEG_mocap_in.chanlocs.tracked_point}');
all_rigidbody = unique(all_tracked_points);

EEG_mocap_out = EEG_mocap_in;
EEG_mocap_out.nbchan = 0;
EEG_mocap_out.chanlocs = [];
EEG_mocap_out.data = [];

for i_rigidbody = 1:length(all_rigidbody)
    
    i_rigidbody
    
    idx_this_rb = find(contains(all_tracked_points,all_rigidbody{i_rigidbody}));
    assert(length(idx_this_rb)==6||length(idx_this_rb)==7,'Rigidbody does not contain 6 or 7 channels.')
    EEG_single_rbs = pop_select( EEG_mocap_in, 'channel',idx_this_rb);
    
    EEG_single_rbs = bemobil_process_rigidbody_data(EEG_single_rbs,lowpass_mocap,lowpass_after_derivative);
    
    EEG_mocap_out.nbchan = EEG_mocap_out.nbchan + EEG_single_rbs.nbchan;
    EEG_mocap_out.chanlocs = [EEG_mocap_out.chanlocs EEG_single_rbs.chanlocs];
    EEG_mocap_out.data = [EEG_mocap_out.data; EEG_single_rbs.data];
end

EEG_mocap_out


% merge rigidbody datasets
EEG_mocap_out = EEG_mocap_in;

EEG_mocap_out.nbchan = EEG_single_rbs(1).nbchan*length(all_rigidbody);
EEG_mocap_out.chanlocs = [EEG_mocap_out.chanlocs EEG_vel.chanlocs EEG_acc.chanlocs];

EEG_mocap_out.data = [EEG_mocap_out.data; EEG_vel.data; EEG_acc.data];
