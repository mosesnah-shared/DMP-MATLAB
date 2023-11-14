function quat_new = LogQuat( quat )
% ===========================================================================
% Descriptions
% ------------
%    Logarithmic Map of unit quaternion
% 
% Parameters
% ----------
%   (1) quat: 1x4 (or 4x1) unit quaternion vector
% 
% Returns
% -------
%   (1) quat_new: conjugated version of quat
%
% ===========================================================================

assert( is_unit_quat( quat ) );

% Create a new quaternion
quat_new = zeros( 'like', quat );

qw = quat_real( quat );
qv = quat_imag( quat );
    
% If qv is zero, then simply zero vector 
% nothing to change
if norm( qv ) >= 1e-6;
    quat_new( 2:4 ) = acos( qw ) * qv/norm( qv );
end

end

