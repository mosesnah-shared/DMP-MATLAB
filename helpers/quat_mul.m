function quat_new = quat_mul( quat1, quat2 )
% ===========================================================================
% Descriptions
% ------------
%    Quaternion multiplication
% 
% Parameters
% ----------
%   (1) quat1: 1x4 (or 4x1) quaternion vector
%   (2) quat2: 1x4 (or 4x1) quaternion vector
% 
% Returns
% -------
%   (1) quat_new = quat1 x quat2
%
% ===========================================================================

assert( ( isrow( quat1 ) && isrow( quat2 ) ) || ( iscolumn( quat1 ) && iscolumn( quat2 ) ) );
assert( ( length( quat1 ) == 4 ) && ( length( quat2 ) == 4 ) )

% For this, return the quaternion as column vector
quat_new = zeros( 4, 1 );

quat_new( 1   ) =  quat1( 1 ) * quat2( 1 ) - dot( quat1( 2 :4 ), quat2( 2 :4 ) );
quat_new( 2:4 ) =  quat1( 1 ) * quat2( 2 : 4 ) + quat2( 1 ) * quat1( 2 : 4 ) + cross( quat1( 2 : 4 ), quat2( 2 : 4 ) );

end

