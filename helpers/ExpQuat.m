function quat_new = ExpQuat( quat )
% ===========================================================================
% Descriptions
% ------------
%    Exponential Map of a quaternion
% 
% Parameters
% ----------
%   (1) quat: 1x4 (or 4x1) quaternion vector
% 
% Returns
% -------
%   (1) quat_new: conjugated version of quat
%
% ===========================================================================

qw = quat_real( quat );
qv = quat_imag( quat );

tmp = norm( qv );

% Create a new quaternion
quat_new = zeros( size( quat ), 'like', quat );

if tmp <= 1e-9
    quat_new( 1   ) = 1;
else
    quat_new( 1   ) = exp( qw ) * cos( tmp );
    quat_new( 2:4 ) = exp( qw ) * sin( tmp )*qv/tmp;
end

end

