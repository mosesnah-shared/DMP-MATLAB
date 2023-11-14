function is_check = is_SO3( R )
% ===========================================================================
% Descriptions
% ------------
%    Check whether the matrix is an element of SO(3)
% 
% Parameters
% ----------
%   (1) mat: 3x3 matrix
% 
% Returns
% -------
%   (1) is_check: boolean value of true/false
%
% ===========================================================================
    
% Matrix should be a 3x3 matrix
assert( all( size( R ) == [ 3, 3 ] ) );

% Threshold value.
thres = 1e-8;

%               RR' = R'R = eye( 3 )                    det( R ) = +1 
is_check = ( norm( R * R' - eye( 3 ) ) <= thres && abs( det( R ) - 1 ) <= thres );

end