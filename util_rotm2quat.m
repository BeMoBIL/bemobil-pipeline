function q = util_rotm2quat( R )
%ROTM2QUAT Convert rotation matrix to quaternion
%   Q = ROTM2QUAT(R) converts a 3D rotation matrix, R, into the corresponding
%   unit quaternion representation, Q. The input, R, is an 3-by-3-by-N matrix
%   containing N orthonormal rotation matrices.
%   The output, Q, is an 4-by-N matrix containing N quaternions. Each
%   quaternion is of the form q = [w x y z], with a scalar number as
%   the first value. Each element of Q must be a real number.
%
%   Example:
%      % Convert a rotation matrix to a quaternion
%      R = [0 0 1; 0 1 0; -1 0 0];
%      q = rotm2quat(R)
%
%   See also util_quat2rotm
%
% Author: Thomas Jespersen, https://github.com/mindThomas/MATLAB-tools/tree/master/Quaternion


% Compute initial quaternion values
qw = sqrt(R(1,1,:)+R(2,2,:)+R(3,3,:)+1)/2.0;
qx = R(3,2,:) - R(2,3,:);
qy = R(1,3,:) - R(3,1,:);
qz = R(2,1,:) - R(1,2,:);

% Compute set of vector elements for general case
qx1 = R(3,1,:) + R(1,3,:);
qy1 = R(3,2,:) + R(2,3,:);
qz1 = R(3,3,:) - R(1,1,:) - R(2,2,:) + 1;
add = (qz >= 0);

% Handle the first degenerate case
case3 = R(2,2,:) >= R(3,3,:);
if any(case3)
    qx1(1,1,case3) = R(2,1,case3) + R(1,2,case3);
    qy1(1,1,case3) = R(2,2,case3) - R(1,1,case3) - R(3,3,case3) + 1;
    qz1(1,1,case3) = R(3,2,case3) + R(2,3,case3);
    add(1,1,case3) = (qy(1,1,case3) >= 0);
end

% Handle the second degenerate case
case2 = (R(1,1,:) >= R(2,2,:)) & (R(1,1,:) >= R(3,3,:));
if any(case2)
    qx1(1,1,case2) = R(1,1,case2) - R(2,2,case2) - R(3,3,case2) + 1;
    qy1(1,1,case2) = R(2,1,case2) + R(1,2,case2);
    qz1(1,1,case2) = R(3,1,case2) + R(1,3,case2);
    add(1,1,case2) = (qx(1,1,case2) >= 0);
end

% Subtract general from initial
qx2 = qx - qx1;
qy2 = qy - qy1;
qz2 = qz - qz1;

% Add in the case that adding is necessary
if any(add)
    qx2(add) = qx(add) + qx1(add);
    qy2(add) = qy(add) + qy1(add);
    qz2(add) = qz(add) + qz1(add);
end

% Compute the norm of the vector components
qnorm = sqrt(bsxfun(@times,qx2,qx2)+bsxfun(@times,qy2,qy2)+bsxfun(@times,qz2,qz2));

zeroNorm = ~qnorm;

% Create a scaling factor for the vector components
if (qw<=1)
qscale = sqrt(1-qw.^2)./qnorm;
else
    qscale=0;
end
% Construct the full quaternion
q = [qw; qscale.*qx2; qscale.*qy2; qscale.*qz2];

% Handle the case of the zero quaternion
if any(zeroNorm)
    q(:,:,~qnorm) = [1; 0; 0; 0];
end

% Shape the output as columns of quaternions
q = reshape(q,[4,numel(q)/4]);

end