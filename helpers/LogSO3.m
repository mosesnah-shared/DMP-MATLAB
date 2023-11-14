function so3 = LogSO3( R )
% ===========================================================================
% Descriptions
% ------------
%    Logarithmic Map of SO3
% 
% Parameters
% ----------
%   (1) R: An SO3 matrix
% 
% Returns
% -------
%   (1) so3: the skew-symmetric part
%
% ===========================================================================

assert( is_SO3( R ) );

% First, need to check whether 
assert( abs( trace( R ) + 1 ) >= 1e-5 );
theta = acos( 0.5 * ( trace( R ) - 1 ) );

if theta <= 1e-3
    so3 = zeros( 3, 3 );
else
    so3 = 1/( 2 * sin( theta ) ) * ( R - R' );
end

end