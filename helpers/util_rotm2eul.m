function eul = util_rotm2eul(rotm, sequence)
% Author: Thomas Jespersen, https://github.com/mindThomas/MATLAB-tools/tree/master/utilities
    if ( (size(rotm,1) ~= 3) || (size(rotm,2) ~= 3) )
        error('rotm2eul: %s', WBM.wbmErrorMsg.WRONG_MAT_DIM);
    end

    if ~exist('sequence', 'var')
        % use the default axis sequence ...
        sequence = 'ZYX';
    end
    eul = zeros(3,1);

    %% Compute the Euler angles theta for the x, y and z-axis from a rotation matrix R, in
    %  dependency of the specified axis rotation sequence for the rotation factorization:
    % For further details see:
    %   [1] Geometric Tools Engine, Documentation: <http://www.geometrictools.com/Documentation/EulerAngles.pdf>, pp. 9-10, 16-17.
    %   [2] Computing Euler angles from a rotation matrix, Gregory G. Slabaugh, <http://www.staff.city.ac.uk/~sbbh653/publications/euler.pdf>
    %   [3] Modelling and Control of Robot Manipulators, L. Sciavicco & B. Siciliano, 2nd Edition, Springer, 2008,
    %       pp. 30-33, formulas (2.19), (2.19'), (2.21) and (2.21').
    switch sequence
        case 'ZYX'
            % convention used by (*) and (**).
            % note: the final orientation is the same as in XYZ order about fixed axes ...
            if (rotm(3,1) < 1)
                if (rotm(3,1) > -1) % case 1: if r31 ~= ±1
                    % Solution with positive sign. It limits the range of the values
                    % of theta_y to (-pi/2, pi/2):
                    eul(1,1) = atan2(rotm(2,1), rotm(1,1)); % theta_z
                    eul(2,1) = asin(-rotm(3,1));            % theta_y
                    eul(3,1) = atan2(rotm(3,2), rotm(3,3)); % theta_x
                else % case 2: if r31 = -1
                    % theta_x and theta_z are linked --> Gimbal lock:
                    % There are infinity number of solutions for theta_x - theta_z = atan2(-r23, r22).
                    % To find a solution, set theta_x = 0 by convention.
                    eul(1,1) = -atan2(-rotm(2,3), rotm(2,2));
                    eul(2,1) = pi/2;
                    eul(3,1) = 0;
                end
            else % case 3: if r31 = 1
                % Gimbal lock: There is not a unique solution for
                %   theta_x + theta_z = atan2(-r23, r22), by convention, set theta_x = 0.
                eul(1,1) = atan2(-rotm(2,3), rotm(2,2));
                eul(2,1) = -pi/2;
                eul(3,1) = 0;
            end
        case 'ZYZ'
            % convention used by (*)
            if (rotm(3,3) < 1)
                if (rotm(3,3) > -1)
                    % Solution with positive sign, i.e. theta_y is in the range (0, pi):
                    eul(1,1) = atan2(rotm(2,3),  rotm(1,3)); % theta_z1
                    eul(2,1) = acos(rotm(3,3));              % theta_y (is equivalent to atan2(sqrt(r13^2 + r23^2), r33) )
                    eul(3,1) = atan2(rotm(3,2), -rotm(3,1)); % theta_z2
                else % if r33 = -1:
                    % Gimbal lock: infinity number of solutions for
                    %   theta_z2 - theta_z1 = atan2(r21, r22), --> set theta_z2 = 0.
                    eul(1,1) = -atan2(rotm(2,1), rotm(2,2)); % theta_z1
                    eul(2,1) = pi;                           % theta_y
                    eul(3,1) = 0;                            % theta_z2
                end
            else % if r33 = 1:
                % Gimbal lock: infinity number of solutions for
                %    theta_z2 + theta_z1 = atan2(r21, r22), --> set theta_z2 = 0.
                eul(1,1) = atan2(rotm(2,1), rotm(2,2)); % theta_z1
                eul(2,1) = 0;                           % theta_y
                eul(3,1) = 0;                           % theta_z2
            end
        % case 'ZYZ-'
        %     % convention used by (**)
        %     if (rotm(3,3) < 1)
        %         if (rotm(3,3) > -1)
        %             % Variant with negative sign. This is a derived solution
        %             % which produces the same effects as the solution above.
        %             % It limits the values of theta_y in the range of (-pi,0):
        %             eul(1,1) = atan2(-rotm(2,3), -rotm(1,3)); % theta_z1
        %             eul(2,1) = -acos(rotm(3,3));              % theta_y (is equivalent to atan2(-sqrt(r13^2 + r23^2), r33) )
        %             eul(3,1) = atan2(-rotm(3,2),  rotm(3,1)); % theta_z2
        %         else % if r33 = -1:
        %             % Gimbal lock: infinity number of solutions for
        %             %   theta_z2 - theta_z1 = atan2(-r12, -r11), --> set theta_z2 = 0.
        %             eul(1,1) = -atan2(-rotm(1,2), -rotm(1,1)); % theta_z1  (correct ???)
        %             eul(2,1) = -pi;                            % theta_y
        %             eul(3,1) = 0;                              % theta_z2
        %         end
        %     else % if r33 = 1:
        %         % Gimbal lock: infinity number of solutions for
        %         %    theta_z2 + theta_z1 = atan2(-r12, -r11), --> set theta_z2 = 0.
        %         eul(1,1) = atan2(-rotm(1,2), -rotm(1,1)); % theta_z1  (correct ???)
        %         eul(2,1) = 0;                             % theta_y
        %         eul(3,1) = 0;                             % theta_z2
        %     end
        otherwise
            error('rotm2eul: %s', WBM.wbmErrorMsg.UNKNOWN_AXIS_SEQ);
    end
end
% (*)  ... The Geometric Tools Engine (http://www.geometrictools.com),
% (**) ... The Robotics System Toolbox for Matlab (http://mathworks.com/help/robotics/index.html).