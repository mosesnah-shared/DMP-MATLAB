function qw = quat_real( quat )
% ===========================================================================
% Descriptions
% ------------
%    Get the imaginary vector part of a quaternion
% 
% Parameters
% ----------
%   (1) quat: 1x4 (or 4x1) quaternion vector
% 
% Returns
% -------
%   (1) qw: Scalar part of the quaternion
%
% ===========================================================================

% The quaternion must either be a column or row vector
assert( isrow( quat ) || iscolumn( quat ) );
assert( length( quat ) == 4 );

qw = quat( 1 );

end

