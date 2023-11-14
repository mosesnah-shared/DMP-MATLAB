function vec = so3_to_R3( mat )
% ===========================================================================
% Descriptions
% ------------
%    Change from so3 matrix to R3 vector
% 
% Parameters
% ----------
%   (1) mat: 3x3 skew-symmetric matrix
% 
% Returns
% -------
%   (1) vec: R3 column vector
%
% ===========================================================================

% Check whether matrix is so3
assert( is_skewsym( mat ) )

% Return a column vector 
vec = [ -mat( 2,3 ); mat( 1,3 ); -mat( 1,2 ) ];

end