% Transforms rigid body dataset unit quaternion angle values (4 dimensions) into eul angles (3 dimensions) to make them
% human-interpretable. Quat values will be taken out, and eul angles will be assigned as new channels.
%
% Input:
%       EEG dataset containing motion channels with quaternion data (must have exactly 4 channels with labels
%       containing "quat_x/y/z/w")
%
% Output:
%       EEG dataset containing motion channels with eul data (will have no quaternion channels but instead
%       eul angle channels with the replaced names of "eul_x/y/z")


function EEG_motion_out = bemobil_motion_quat2eul(EEG_motion_in)


for channel = 1:EEG_motion_in.nbchan
    
    % checking for already present euls
    if any(~cellfun(@isempty,strfind(lower({EEG_motion_in.chanlocs.labels}),'eul')))
        error('Dataset already contains eul data.')
    end
    
end

assert(EEG_motion_in.nbchan == 7,'Rigidbody quaternion dataset needs to have exactly 7 channels: XYZ position and XYZW quaternions.')
% no euls are present, therefore each RB has 7 channels:
% XYZABCD, from which ABCD are the quaternion values

data = EEG_motion_in.data';
newData = zeros(size(data,1),6);
newLabel = cell(6,1);
% the new eul data has 1 channel less than the quaternions

% fill the new data set and its label with all initial position data
channel_labels = {EEG_motion_in.chanlocs.labels};
non_quat_indices = cellfun(@isempty,strfind(lower(channel_labels),'quat'));

newLabel(1:sum(non_quat_indices)) = channel_labels(non_quat_indices);
newData(:,1:sum(non_quat_indices)) = data(:,non_quat_indices);

% now fill with eul data
quat_indices = ~cellfun(@isempty,strfind(lower(channel_labels),'quat'));
assert(sum(quat_indices)==4,'There must be exactly 4 quaternion channels containing the label ''quat_<x,y,z,w>''!')


% find correct channelnumber for the quaternion values of
% this RB
quaternionX = ~cellfun(@isempty,strfind(lower(channel_labels),'quat_x'));
assert(sum(quaternionX)==1,'There must be exactly 1 quaternion channel containing the label ''quat_x''!')
quaternionY = ~cellfun(@isempty,strfind(lower(channel_labels),'quat_y'));
assert(sum(quaternionY)==1,'There must be exactly 1 quaternion channel containing the label ''quat_y''!')
quaternionZ = ~cellfun(@isempty,strfind(lower(channel_labels),'quat_z'));
assert(sum(quaternionZ)==1,'There must be exactly 1 quaternion channel containing the label ''quat_z''!')
quaternionW = ~cellfun(@isempty,strfind(lower(channel_labels),'quat_w'));
assert(sum(quaternionW)==1,'There must be exactly 1 quaternion channel containing the label ''quat_w''!')

% take the values
x = data(:,quaternionX);
y = data(:,quaternionY);
z = data(:,quaternionZ);
w = data(:,quaternionW);

% check if values are [-1 1] - could have been messed up by
% interpolating or filtering, would be undefined then.
w(w>=1) = 0.99999;
w(w<=-1) = -0.99999;
x(x>=1) = 0.99999;
x(x<=-1) = -0.99999;
y(y>=1) = 0.99999;
y(y<=-1) = -0.99999;
z(z>=1) = 0.99999;
z(z<=-1) = -0.99999;

% transform to ZYX sequence eul angles
disp('Transforming Quaternions to eul angles using a Body-ZYX sequence!')
euls = util_quat2eul([w x y z]');

% convert from radian to degree
euls = euls'*180/pi;

% fill new data set and labels
newData(:,4:6) = euls;

% take the original prefix before 'quat_x' as a prefix for all new eul channels
newLabel{4} = strcat(channel_labels{quaternionX}(1:strfind(lower(channel_labels{quaternionX}),'quat_x')-1),'eul_x');
newLabel{5} = strcat(channel_labels{quaternionX}(1:strfind(lower(channel_labels{quaternionX}),'quat_x')-1),'eul_y');
newLabel{6} = strcat(channel_labels{quaternionX}(1:strfind(lower(channel_labels{quaternionX}),'quat_x')-1),'eul_z');

% make new set
EEG_motion_out = EEG_motion_in;
EEG_motion_out.nbchan = 6;
EEG_motion_out.data = newData';
EEG_motion_out.chanlocs(end) = [];
EEG_motion_out.chanlocs(1).labels = newLabel{1};
EEG_motion_out.chanlocs(2).labels = newLabel{2};
EEG_motion_out.chanlocs(3).labels = newLabel{3};
EEG_motion_out.chanlocs(4).labels = newLabel{4};
EEG_motion_out.chanlocs(5).labels = newLabel{5};
EEG_motion_out.chanlocs(6).labels = newLabel{6};

EEG_motion_out.etc.quat2eul_sequence = 'Body-ZYX';
