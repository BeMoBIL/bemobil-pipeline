function q = util_eul2quat( eul, varargin )
%EUL2QUAT Convert Euler angles to quaternion
%   Q = EUL2QUAT(EUL) converts a given set of 3D Euler angles, EUL, into
%   the corresponding unit quaternion, Q. EUL is an 3-by-N matrix of Euler 
%   rotation angles.
%   The output, Q, is an 4-by-N matrix containing N quaternions. Each
%   quaternion is of the form q = [w x y z]', with a scalar number as
%   the first value.
%
%   Q = EUL2QUAT(EUL, SEQ) converts a set of 3D Euler angles into a unit
%   quaternion. The Euler angles are specified by the axis rotation
%   sequence, SEQ.
%
%   The default rotation sequence is 'ZYX', where the order of rotation
%   angles is Z Axis Rotation, Y Axis Rotation, and X Axis Rotation.
%
%   The following rotation sequences, SEQ, are supported: 'ZYX' and 'ZYZ'.
%
%   Example:
%      % Calculate the quaternion for a set of Euler angles
%      % By default, the ZYX axis order will be used.
%      angles = [0 pi/2 0];
%      q = eul2quat(angles)
%
%      % Calculate the quaternion based on a ZYZ rotation
%      qzyz = eul2quat(angles, 'ZYZ')
%
%   See also util_quat2eul
%
% Author: Thomas Jespersen, https://github.com/mindThomas/MATLAB-tools/tree/master/Quaternion

if (length(varargin) == 1)
    seq = varargin{1};
else
    seq = 'ZYX';
end

if ( (size(eul,1) == 1) && (size(eul,2) == 3) )
    eul = eul';
    transpose = true;
else
    transpose = false;
end

% Compute sines and cosines of half angles
c = cos(eul/2);
s = sin(eul/2);

% The parsed sequence will be in all upper-case letters and validated
switch seq
    case 'ZYX'
        % Construct quaternion
        q = [c(1,:).*c(2,:).*c(3,:)+s(1,:).*s(2,:).*s(3,:); ...
            c(1,:).*c(2,:).*s(3,:)-s(1,:).*s(2,:).*c(3,:); ...
            c(1,:).*s(2,:).*c(3,:)+s(1,:).*c(2,:).*s(3,:); ...
            s(1,:).*c(2,:).*c(3,:)-c(1,:).*s(2,:).*s(3,:)];
        
    case 'ZYZ'
        % Construct quaternion
        q = [c(1,:).*c(2,:).*c(3,:)-s(1,:).*c(2,:).*s(3,:); ...
            c(1,:).*s(2,:).*s(3,:)-s(1,:).*s(2,:).*c(3,:); ...
            c(1,:).*s(2,:).*c(3,:)+s(1,:).*s(2,:).*s(3,:); ...
            s(1,:).*c(2,:).*c(3,:)+c(1,:).*c(2,:).*s(3,:)];
end

if (transpose)
    q = q';
end

end