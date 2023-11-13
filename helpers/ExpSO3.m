function R = ExpSO3( so3 )

assert( is_skewsym( so3 ) );
theta = norm( so3_to_R3( so3 ) );
R = eye( 3 ) + sin( theta )/theta * so3 + ( 1 - cos( theta ) )/theta^2 * so3^2;

end