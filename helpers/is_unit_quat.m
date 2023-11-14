function is_check = is_unit_quat( quat )
% ===========================================================================
% Descriptions
% ------------
%    Check whether the the 4D vector is a unit quaternion
% 
% Parameters
% ----------
%   (1) quat: 1x4 (or 4x1) unit quaternion
% 
% Returns
% -------
%   (1) is_check: boolean value of true/false
%
% ===========================================================================
    
% The quaternion vector must be a row or column vector
assert( isrow( quat ) || iscolumn( quat )  )

% Set the threshold 
thres = 1e-8;

% Get the no
is_check = ( abs( norm( quat ) - 1 ) <= thres );

end