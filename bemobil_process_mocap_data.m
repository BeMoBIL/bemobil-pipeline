% Process a dataset containing rigid body data streams (which contain 6 DOF data of position and orientation).
% Orientation can be either euler angles or unit quaternions. Orientation channels need to have the suffixes 'eul_x/y/z'
% or 'quat_x/y/z/w'. The dataset must contain information about the rigid bodies in the chanlocs substruct, specifically
% the respective rigid bodies need to be marked by the same "EEG.chanlocs.tracked_location" entry. Processing contains
% optional lowpass filtering, transformation to euler angles, and taking the first two derivatives (velocity and
% acceleration).
%
% Input:
%       EEG_mocap_in                - EEGLAB dataset containing mocap rigid body data (6 DOF)
%       lowpass_mocap               - (OPTIONAL) lowpass passband edge for filtering the raw mocap data. If empty, no 
%                                       filter is applied. 7Hz is recommended.
%       lowpass_after_derivative    - (OPTIONAL) lowpass passband edge for filtering during after time derivative. If 
%                                       empty, no filter is applied. Can help if the derivatives are too noisy, but not
%                                       specifically recommended.
%
% Output:
%       EEG_mocap_out               - EEGLAB dataset containing processed rigid body data with 6 DOF in Euler angles,
%                                       plus the first two derivatives (18 channels altogether)

function EEG_mocap_out = bemobil_process_mocap_data(EEG_mocap_in,lowpass_mocap,lowpass_after_derivative)

if ~exist('lowpass_mocap','var') || isempty(lowpass_mocap)
    lowpass_mocap = [];
end
if ~exist('lowpass_after_derivative','var') || isempty(lowpass_after_derivative)
    lowpass_after_derivative = [];
end

all_tracked_points = lower({EEG_mocap_in.chanlocs.tracked_point}');
all_rigidbody = unique(all_tracked_points);

EEG_mocap_out = EEG_mocap_in;
EEG_mocap_out.nbchan = 0;
EEG_mocap_out.chanlocs = [];
EEG_mocap_out.data = [];

for i_rigidbody = 1:length(all_rigidbody)
    
    disp('-----------------------------------------------------------------')
    disp([num2str(i_rigidbody) '/' num2str(length(all_rigidbody)) ': Processing ' all_rigidbody{i_rigidbody} ' rigid body data!'])
    disp('-----------------------------------------------------------------')
    
    idx_this_rb = find(contains(all_tracked_points,all_rigidbody{i_rigidbody}));
    assert(length(idx_this_rb)==6||length(idx_this_rb)==7,'Rigidbody does not contain 6 (Euler) or 7 (quaternions) channels!')
    EEG_single_rbs = pop_select( EEG_mocap_in, 'channel',idx_this_rb);
    
    EEG_single_rbs = bemobil_process_rigidbody_data(EEG_single_rbs,lowpass_mocap,lowpass_after_derivative);
    
    EEG_mocap_out.nbchan = EEG_mocap_out.nbchan + EEG_single_rbs.nbchan;
    EEG_mocap_out.chanlocs = [EEG_mocap_out.chanlocs EEG_single_rbs.chanlocs];
    EEG_mocap_out.data = [EEG_mocap_out.data; EEG_single_rbs.data];
end

