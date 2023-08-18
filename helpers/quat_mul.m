function quat_new = quat_mul( quat1, quat2 )

assert( all( size( quat1 ) == [ 1,4 ] ) || all( size( quat1 ) == [ 4,1 ] ) );
assert( all( size( quat2 ) == [ 1,4 ] ) || all( size( quat2 ) == [ 4,1 ] ) );

            
if all( [ 1,4 ] == size( quat1 )  )
    quat1 = quat1'; 
end

if all( [ 1,4 ] == size( quat2 )  )
    quat2 = quat2'; 
end

quat_new = zeros( 4, 1 );

quat_new( 1   ) =  quat1( 1 ) * quat2( 1 ) - dot( quat1( 2 :4 ), quat2( 2 :4 ) );
quat_new( 2:4 ) =  quat1( 1 ) * quat2( 2 : 4 ) + quat2( 1 ) * quat1( 2 : 4 ) + cross( quat1( 2 : 4 ), quat2( 2 : 4 ) );

end

