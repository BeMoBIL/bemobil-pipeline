function EEG_mocap_out = bemobil_mocap_eul2quat(EEG_mocap_in)
% Transforms Euler angle values (3 dimensions) into unit quaternions (4 dimensions).
% Euler angles will be taken out, and Quat values will be assigned as new channels
%
% Input arguments:
%       EEG dataset containing mocap channels with euler data (must have 3 euler channels with labels containing
%       "euler_x/y/z")
%
%
% Output argument:
%       EEG dataset containing mocap channels with quaternion data (will have 3 channels with labels containing
%       "quat_x/y/z/w")
%
% Usage:


for channel = 1:EEG_mocap_in.nbchan
    
    % checking for already present eulers
    if any(~cellfun(@isempty,strfind(lower({EEG_mocap_in.chanlocs.labels}),'quat')))
        error('Dataset already contains quaternion data.')
    end
    
end

assert(EEG_mocap_in.nbchan == 6,'Rigidbody euler angle dataset needs to have exactly 6 channels: XYZ position and XYZ orientation.')


data = EEG_mocap_in.data';
newData = zeros(size(data,1),7);
newLabel = cell(7,1);
% the new quaternion data has 1 channel more than the quaternions

% fill the new data set and its label with all initial position data
channel_labels = {EEG_mocap_in.chanlocs.labels};
non_eul_indices = cellfun(@isempty,strfind(lower(channel_labels),'euler'));

newLabel(1:sum(non_eul_indices)) = channel_labels(non_eul_indices);
newData(:,1:sum(non_eul_indices)) = data(:,non_eul_indices);

% now fill with euler data
quat_indices = ~cellfun(@isempty,strfind(lower(channel_labels),'euler'));
assert(sum(quat_indices)==3,'There must be exactly 3 euler angle channels containing the label ''euler_<x,y,z,w>''!')


% find correct channelnumber for the quaternion values of
% this RB
quaternionX = ~cellfun(@isempty,regexp(lower(channel_labels),'euler_x'));
assert(sum(quaternionX)==1,'There must be exactly 1 quaternion channel containing the label ''euler_x''!')
quaternionY = ~cellfun(@isempty,regexp(lower(channel_labels),'euler_y'));
assert(sum(quaternionY)==1,'There must be exactly 1 quaternion channel containing the label ''euler_y''!')
quaternionZ = ~cellfun(@isempty,regexp(lower(channel_labels),'euler_z'));
assert(sum(quaternionZ)==1,'There must be exactly 1 quaternion channel containing the label ''euler_z''!')

% take the values
x = data(:,quaternionX);
y = data(:,quaternionY);
z = data(:,quaternionZ);

eulers = [x y z];


if any([range(eulers)>2+pi])
    % convert to radians if is in degrees
    eulers = eulers/(180/pi);
end

% transform to ZYX sequence euler angles
disp('Transforming Quaternions to Euler angles using a Body-ZYX sequence!')
quats = util_eul2quat(eulers');

% fill new data set and labels
newData(:,4:7) = quats';

% take the original prefix before 'quat_x' as a prefix for all new euler channels
newLabel{4} = strcat(channel_labels{4}(1:strfind(lower(channel_labels{4}),'euler_x')-1),'quat_w');
newLabel{5} = strcat(channel_labels{4}(1:strfind(lower(channel_labels{4}),'euler_x')-1),'quat_x');
newLabel{6} = strcat(channel_labels{4}(1:strfind(lower(channel_labels{4}),'euler_x')-1),'quat_y');
newLabel{7} = strcat(channel_labels{4}(1:strfind(lower(channel_labels{4}),'euler_x')-1),'quat_z');

% make new set
EEG_mocap_out = EEG_mocap_in;
EEG_mocap_out.nbchan = 7;
EEG_mocap_out.data = newData';
EEG_mocap_out.chanlocs(end) = [];
EEG_mocap_out.chanlocs(4).labels = newLabel{4};
EEG_mocap_out.chanlocs(5).labels = newLabel{5};
EEG_mocap_out.chanlocs(6).labels = newLabel{6};
EEG_mocap_out.chanlocs(7).labels = newLabel{7};

EEG_mocap_out.etc.eul2quat_sequence = 'Body-ZYX';
