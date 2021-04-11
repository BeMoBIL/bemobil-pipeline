function EEG_mocap_out = bemobil_process_rigidbody_data(EEG_mocap_in,lowpass_mocap,lowpass_after_derivative)
% Process a single rigid body data stream which contains 6 DOF data of position and orientation. Orientation can be
% either eul angles. Channels need to have the suffixes 'eul_x/y/z'), or unit quaternions (channels need to have the
% suffixes 'quat_x/y/z/w'. Processing contains lowpass filtering, transformation to eul angles, and taking the first
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
%       EEG_mocap_out   - EEGLAB dataset containing processed mocap channels with 6 DOF in eul angles, plus the first
%       two derivatives (18 channels altogether)

% make sure to use double precision
try
    pop_editoptions( 'option_single', 0);
catch
    warning('Tried editing EEGLAB options to use double precision but failed! Maybe the file was in use by another MATLAB instance?')
end

% make sure all euler values are not exceeding pi due to some weirdness in the upsampling

eul_indices = ~cellfun(@isempty,strfind(lower({EEG_mocap_in.chanlocs.labels}),'eul'));
euldata = EEG_mocap_in.data(eul_indices,:);
euldata(euldata>pi)=pi;
euldata(euldata<-pi)=-pi;
EEG_mocap_in.data(eul_indices,:) = euldata;

% empty EEGLAB entries
ALLEEG = []; EEG=[]; CURRENTSET=[];

% make sure data is in quaternion angles, so we can filter the data without jumps from 0 to 360 degrees
try
    EEG = bemobil_mocap_eul2quat(EEG_mocap_in);
catch ME
    warning(ME.message)
    if strcmp(ME.message,'You can only unflip Quaternions, this dataset contains eul angles, try it with the original data set.')
        disp('Attempt to transform from eul to quaternion angles failed. Assuming original dataset is in quaternion angles.')
        EEG = EEG_mocap_in;
    else
        error('Unexpected error (see warning above).')
    end
end

% unflip quaternion angles if they are flipped (mathematically there are two way to represent the same angle in
% quaternion units, which is totally fine, but if they are flipping we cannot filter the data as there would be
% distortions and ringing artifacts
EEG = bemobil_mocap_unflipSigns(EEG);

% lowpass mocap data
if exist('lowpass_mocap','var') && ~isempty(lowpass_mocap)
    [ ALLEEG EEG CURRENTSET ] = bemobil_filter(ALLEEG, EEG, CURRENTSET, [],lowpass_mocap,[], []);
end

% transform back to eul angles to make it interpretable. Using the "Body-ZYX" scheme.
EEG = bemobil_mocap_quat2eul(EEG);

% create new derivatives datasets (filtering after taking a derivative is recommended but optional, as just using the
% diff increases the noise level). Derivatives account for jumps in eul orientation from 0 to 360 degrees and subtract
% accordingly.
EEG_vel = bemobil_mocap_timeDerivative(EEG);
if exist('lowpass_after_derivative','var') && ~isempty(lowpass_after_derivative)
    [ ALLEEG EEG_vel CURRENTSET ] = bemobil_filter(ALLEEG, EEG_vel, CURRENTSET, [],lowpass_after_derivative,[], []);
end

EEG_acc = bemobil_mocap_timeDerivative(EEG_vel);
if exist('lowpass_after_derivative','var') && ~isempty(lowpass_after_derivative)
    [ ALLEEG EEG_acc CURRENTSET ] = bemobil_filter(ALLEEG, EEG_acc, CURRENTSET, [],lowpass_after_derivative,[], []);
end

% merge the three mocap datasets
EEG_mocap_out = EEG;

EEG_mocap_out.nbchan = 18;
EEG_mocap_out.chanlocs = [EEG_mocap_out.chanlocs EEG_vel.chanlocs EEG_acc.chanlocs];

EEG_mocap_out.data = [EEG_mocap_out.data; EEG_vel.data; EEG_acc.data];
