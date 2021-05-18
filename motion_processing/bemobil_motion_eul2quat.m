function EEG_motion_out = bemobil_motion_eul2quat(EEG_motion_in)
% Transforms eul angle values (3 dimensions) into unit quaternions (4 dimensions).
% eul angles will be taken out, and Quat values will be assigned as new channels
%
% Input arguments:
%       EEG dataset containing motion channels with eul data (must have 3 eul channels with labels containing
%       "eul_x/y/z")
%
%
% Output argument:
%       EEG dataset containing motion channels with quaternion data (will have 3 channels with labels containing
%       "quat_x/y/z/w")
%
% Usage:


for channel = 1:EEG_motion_in.nbchan
    
    % checking for already present euls
    if any(~cellfun(@isempty,strfind(lower({EEG_motion_in.chanlocs.labels}),'quat')))
        error('Dataset already contains quaternion data.')
    end
    
end

assert(EEG_motion_in.nbchan == 6,'Rigidbody eul angle dataset needs to have exactly 6 channels: XYZ position and XYZ orientation.')


data = EEG_motion_in.data';
newData = zeros(size(data,1),7);
newLabel = cell(7,1);
% the new quaternion data has 1 channel more than the quaternions

% fill the new data set and its label with all initial position data
channel_labels = {EEG_motion_in.chanlocs.labels};
non_eul_indices = cellfun(@isempty,strfind(lower(channel_labels),'eul'));

newLabel(1:sum(non_eul_indices)) = channel_labels(non_eul_indices);
newData(:,1:sum(non_eul_indices)) = data(:,non_eul_indices);

% now fill with eul data
eul_indices = ~cellfun(@isempty,strfind(lower(channel_labels),'eul'));
assert(sum(eul_indices)==3,'There must be exactly 3 eul angle channels containing the label ''eul_<x,y,z,w>''!')


% find correct channelnumber for the quaternion values of
% this RB
eulX = ~cellfun(@isempty,regexp(lower(channel_labels),'eul_x'));
assert(sum(eulX)==1,'There must be exactly 1 quaternion channel containing the label ''eul_x''!')
eulY = ~cellfun(@isempty,regexp(lower(channel_labels),'eul_y'));
assert(sum(eulY)==1,'There must be exactly 1 quaternion channel containing the label ''eul_y''!')
eulZ = ~cellfun(@isempty,regexp(lower(channel_labels),'eul_z'));
assert(sum(eulZ)==1,'There must be exactly 1 quaternion channel containing the label ''eul_z''!')

% take the values
x = data(:,eulX);
y = data(:,eulY);
z = data(:,eulZ);

euls = [x y z];


if any([range(euls)>2+pi])
    % convert to radians if is in degrees
    euls = euls/(180/pi);
end

% transform to ZYX sequence eul angles
disp('Transforming Euler to Quaternion angles using a Body-ZYX sequence!')
quats = util_eul2quat(euls');

% fill new data set and labels
newData(:,4:7) = quats';

% take the original prefix before 'quat_x' as a prefix for all new eul channels
newLabel{4} = strcat(channel_labels{eulX}(1:strfind(lower(channel_labels{eulX}),'eul_x')-1),'quat_w');
newLabel{5} = strcat(channel_labels{eulX}(1:strfind(lower(channel_labels{eulX}),'eul_x')-1),'quat_x');
newLabel{6} = strcat(channel_labels{eulX}(1:strfind(lower(channel_labels{eulX}),'eul_x')-1),'quat_y');
newLabel{7} = strcat(channel_labels{eulX}(1:strfind(lower(channel_labels{eulX}),'eul_x')-1),'quat_z');

% make new set
EEG_motion_out = EEG_motion_in;
EEG_motion_out.nbchan = 7;
EEG_motion_out.data = newData';
EEG_motion_out.chanlocs(end) = [];
EEG_motion_out.chanlocs(1).labels = newLabel{1};
EEG_motion_out.chanlocs(2).labels = newLabel{2};
EEG_motion_out.chanlocs(3).labels = newLabel{3};
EEG_motion_out.chanlocs(4).labels = newLabel{4};
EEG_motion_out.chanlocs(5).labels = newLabel{5};
EEG_motion_out.chanlocs(6).labels = newLabel{6};
EEG_motion_out.chanlocs(7).labels = newLabel{7};

EEG_motion_out.etc.eul2quat_sequence = 'Body-ZYX';
