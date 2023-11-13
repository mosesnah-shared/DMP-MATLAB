function so3 = LogSO3( R )

assert( is_SO3( R ) );

% Get the theta value 
% If trace value 
assert( abs( trace( R ) + 1 ) >= 1e-5 );
theta = acos( 0.5 * ( trace( R ) - 1 ) );

if theta <= 1e-3
    so3 = zeros( 3, 3 );
else
    so3 = 1/( 2 * sin( theta ) ) * ( R - R' );
end

end