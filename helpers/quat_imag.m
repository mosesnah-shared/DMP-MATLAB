function vec = quat_imag( quat )
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
%   (1) vec: the 3D vector part of quaternion
%
% ===========================================================================

% The quaternion must either be a column or row vector
assert( isrow( quat ) || iscolumn( quat ) );
assert( length( quat ) == 4 );

vec = quat( 2:4 ); 

end

