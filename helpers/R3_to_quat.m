function quat = R3_to_quat( vec )
% ===========================================================================
% Descriptions
% ------------
%    Change from R3 vector to quaternion vector
% 
% Parameters
% ----------
%   (1) vec: 3x1 (or 1x3) vector
% 
% Returns
% -------
%   (1) so3: the skew-symmetric form of vec
%
% ===========================================================================

% Should either be a column or row vector
assert( iscolumn( vec ) || isrow( vec ) );

% The length should be three
assert( length( vec ) == 3 );

quat = zeros( 1, 4 );
quat( 2:4 ) = vec;
         
end