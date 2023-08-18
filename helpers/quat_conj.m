function quat_new = quat_conj( quat )

assert( all( size( quat ) == [ 1,4 ] ) || all( size( quat ) == [ 4,1 ] ) );

quat_new = zeros( 'like', quat );

quat_new( 1   ) =  quat( 1 );
quat_new( 2:4 ) = -quat( 2:4 );

end

