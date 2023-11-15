function vec = get_quat_error( quat1, quat2 )
% ===========================================================================
% Descriptions
% ------------
%    Get the error 3D vector between two unit quaternions
% 
% Parameters
% ----------
%   (1) quat1: 1x4 (or 4x1) unit quaternion vector
%   (2) quat2: 1x4 (or 4x1) unit quaternion vector
% 
% Returns
% -------
%   (1) vec: 2 * Im( Log( quat1^* x quat2 )  )
%
% ===========================================================================

% The quaternion must either be a column or row vector
assert( is_unit_quat( quat1 ) && is_unit_quat( quat2 ) )

vec = 2 * quat_imag( LogQuat(  quat_mul( quat_conj( quat1 ), quat2 )  )  );

end

