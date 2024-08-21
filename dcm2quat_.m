function q = dcm2quat_( dcm, action, tol )
%  DCM2QUAT Convert direction cosine matrix to quaternion.
%   Q = DCM2QUAT( DCM ) calculates the quaternion, Q, for given direction
%   cosine matrix, DCM. 
% 
%   Q = DCM2QUAT( DCM, ACTION ) uses ACTION for error handling during
%   rotation matrix validation. 
%
%   Q = DCM2QUAT( DCM, ACTION, TOL ) uses TOL as a relative error tolerance
%   for rotation matrix validation. 
%   eps(2).
%
% Inputs:
%   DCM    :3-by-3-by-M matrix containing M orthogonal direction cosine matrices.
%           DCM performs the coordinate transformation from inertial axes to body axes. 
%   ACTION :Scalar string specifying the action to take if the DCM is invalid.
%           Valid directon cosine matricies are orthogonal and proper:
%           transpose(DCM) * DCM == 1 +/- TOL
%           det(DCM) == 1 +/- TOL
%           ACTION must be one of the following:
%             'Error'
%             'Warning'
%             'None'
%           Default is 'None'.
%   TOL    :Scalar numeric specifying the relative tolerance for rotation matrix validation.
%           The tolerance is used for validating the direction cosine matricies.
%           Default value is eps(2).
%
% Outputs:
%   Q :M-by-4 matrix of quaternions vectors. Q has its scalar number as the first
%      column.
%
% Examples:
%   Determine the quaternion from a direction cosine matrix:
%      dcm = [0 1 0; 1 0 0; 0 0 -1];
%      q = dcm2quat(dcm)
%
%   Determine the quaternions from multiple direction cosine matrices:
%      dcm        = [ 0 1 0; 1 0 0; 0 0 -1]; 
%      dcm(:,:,2) = [ 0.4330    0.2500   -0.8660; ...
%                     0.1768    0.9186    0.3536; ...
%                     0.8839   -0.3062    0.3536];
%      q = dcm2quat(dcm)
%
%   Determine the quaternion from a direction cosine matrix with a specified
%   action and tolerance:
%      dcm = [0 1 0; 1 0 0; 0 0 -1];
%      q = dcm2quat(dcm, 'warning', 0.01)
%
%   See also ANGLE2DCM, DCM2ANGLE, QUAT2DCM, QUAT2ANGLE, ANGLE2QUAT.

%   Copyright 2000-2020 The MathWorks, Inc.

for i = size(dcm,3):-1:1

    q(i,4) =  0; 
    
    tr = trace(dcm(:,:,i));

    if (tr > 0)
        sqtrp1 = sqrt( tr + 1.0 );
        
        q(i,1) = 0.5*sqtrp1; 
        q(i,2) = (dcm(2, 3, i) - dcm(3, 2, i))/(2.0*sqtrp1);
        q(i,3) = (dcm(3, 1, i) - dcm(1, 3, i))/(2.0*sqtrp1); 
        q(i,4) = (dcm(1, 2, i) - dcm(2, 1, i))/(2.0*sqtrp1); 
    else
        d = diag(dcm(:,:,i));
        if ((d(2) > d(1)) && (d(2) > d(3)))
            % max value at dcm(2,2,i)
            sqdip1 = sqrt(d(2) - d(1) - d(3) + 1.0 );
            
            q(i,3) = 0.5*sqdip1; 
            
            if ( sqdip1 ~= 0 )
                sqdip1 = 0.5/sqdip1;
            end
            
            q(i,1) = (dcm(3, 1, i) - dcm(1, 3, i))*sqdip1; 
            q(i,2) = (dcm(1, 2, i) + dcm(2, 1, i))*sqdip1; 
            q(i,4) = (dcm(2, 3, i) + dcm(3, 2, i))*sqdip1; 
        elseif (d(3) > d(1))
            % max value at dcm(3,3,i)
            sqdip1 = sqrt(d(3) - d(1) - d(2) + 1.0 );
            
            q(i,4) = 0.5*sqdip1; 
            
            if ( sqdip1 ~= 0 )
                sqdip1 = 0.5/sqdip1;
            end
            
            q(i,1) = (dcm(1, 2, i) - dcm(2, 1, i))*sqdip1;
            q(i,2) = (dcm(3, 1, i) + dcm(1, 3, i))*sqdip1; 
            q(i,3) = (dcm(2, 3, i) + dcm(3, 2, i))*sqdip1; 
        else
            % max value at dcm(1,1,i)
            sqdip1 = sqrt(d(1) - d(2) - d(3) + 1.0 );
            
            q(i,2) = 0.5*sqdip1; 
            
            if ( sqdip1 ~= 0 )
                sqdip1 = 0.5/sqdip1;
            end
            
            q(i,1) = (dcm(2, 3, i) - dcm(3, 2, i))*sqdip1; 
            q(i,3) = (dcm(1, 2, i) + dcm(2, 1, i))*sqdip1; 
            q(i,4) = (dcm(3, 1, i) + dcm(1, 3, i))*sqdip1; 
        end
    end
end

end