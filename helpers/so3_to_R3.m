function vec = so3_to_R3( mat )

assert( all( size( mat ) == [ 3, 3 ] ) )
assert( is_skewsym( mat ) );

vec = [ -mat( 2,3 ); mat( 1,3 ); -mat( 1,2 ) ];

end