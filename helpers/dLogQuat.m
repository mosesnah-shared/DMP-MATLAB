function vec_new = dLogQuat( quat, dquat )
% ===========================================================================
% Descriptions
% ------------
%    Derivative of the Logarithmic Map of unit quaternion
% 
% Parameters
% ----------
%   (1) quat: 1x4 (or 4x1) unit quaternion vector
%
%   (2) dquat: time derivative of the quaternion
% 
% Returns
% -------
%   (1) vec: 1x3 result
%
% ===========================================================================

% Check unit quaternion
assert( is_unit_quat( quat ) );

% dquat should be a column vector for convenience
assert( iscolumn( dquat ) );

% Getting out the parts
eta = quat_real( quat );
eps = quat_imag( quat );

% The Jacobian matrix
a1 = ( -norm( eps ) + acos( eta ) * eta )/norm( eps )^3;
a2 = acos( eta )/norm( eps );
JQ = [ a1 * eps, eye( 3 ) * a2 ];

vec_new = JQ * dquat;

end

