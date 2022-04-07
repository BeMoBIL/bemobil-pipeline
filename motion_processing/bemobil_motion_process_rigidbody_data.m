% Process a single rigid body data stream which contains 6 DOF data of position and orientation. Orientation can be
% either eul angles or unit quaternions. Orientation channels need to have the suffixes 'eul_x/y/z' or 'quat_x/y/z/w'.
% Processing contains lowpass filtering, transformation to eul angles, and taking the first two derivatives (velocity
% and acceleration).
%
% Input:
%       EEG_motion_in                - EEGLAB dataset containing motion channels with 6 DOF
%       lowpass_motion               - (OPTIONAL) lowpass passband edge for filtering the raw motion data. If empty, no 
%                                       filter is applied. 7Hz is recommended.
%       lowpass_after_derivative    - (OPTIONAL) lowpass passband edge for filtering during after time derivative. If 
%                                       empty, no filter is applied. Can help if the derivatives are too noisy, but not
%                                       specifically recommended.
%
% Output:
%       EEG_motion_out               - EEGLAB dataset containing processed motion channels with 6 DOF in eul angles, plus
%                                       the first two derivatives (18 channels altogether)

function EEG_motion_out = bemobil_motion_process_rigidbody_data(EEG_motion_in,lowpass_motion,lowpass_after_derivative)

% make sure to use double precision
try
    pop_editoptions( 'option_single', 0);
catch
    warning('Tried editing EEGLAB options to use double precision but failed! Maybe the file was in use by another MATLAB instance?')
end

% make sure data is in euler angles, so we can clean the data
try
    EEG = bemobil_motion_quat2eul(EEG_motion_in);
catch ME
    if strcmp(ME.message,'Dataset already contains eul data.')
        disp('Assuming original dataset is in euler angles.')
        EEG = EEG_motion_in;
    else
        error('Unexpected error see following warning:')
        warning(ME.message)
    end
end

% clean euler angle data
EEG = bemobil_clean_motion(EEG,20);

% empty EEGLAB entries
ALLEEG = []; CURRENTSET=[];

% unwrap the euler angles to enable filtering
eul_indices = find(~cellfun(@isempty,strfind(lower({EEG.chanlocs.labels}),'eul')));
EEG.data(eul_indices,:) = unwrap(EEG.data(eul_indices,:),[],2);

% lowpass motion data
if exist('lowpass_motion','var') && ~isempty(lowpass_motion)
    [ ALLEEG EEG CURRENTSET ] = bemobil_filter(ALLEEG, EEG, CURRENTSET, [],lowpass_motion,[], []);
end

% wrap the euler angles back to pi
for idx_chan = eul_indices
    EEG.data(idx_chan,:) = wrapToPi(EEG.data(idx_chan,:));
end

% create new derivatives datasets (filtering after taking a derivative is recommended but optional, as just using the
% diff increases the noise level). Derivatives account for jumps in eul orientation from 0 to 360 degrees and subtract
% accordingly.
EEG_vel = bemobil_motion_timeDerivative(EEG);

% filter if set up
if exist('lowpass_after_derivative','var') && ~isempty(lowpass_after_derivative)
    [ ALLEEG EEG_vel CURRENTSET ] = bemobil_filter(ALLEEG, EEG_vel, CURRENTSET, [],lowpass_after_derivative,[], []);
end

EEG_acc = bemobil_motion_timeDerivative(EEG_vel);

% filter if set up
if exist('lowpass_after_derivative','var') && ~isempty(lowpass_after_derivative)
    [ ALLEEG EEG_acc CURRENTSET ] = bemobil_filter(ALLEEG, EEG_acc, CURRENTSET, [],lowpass_after_derivative,[], []);
end

% merge the three motion datasets
EEG_motion_out = EEG;

EEG_motion_out.nbchan = 18;
EEG_motion_out.chanlocs = [EEG_motion_out.chanlocs EEG_vel.chanlocs EEG_acc.chanlocs];

EEG_motion_out.data = [EEG_motion_out.data; EEG_vel.data; EEG_acc.data];
