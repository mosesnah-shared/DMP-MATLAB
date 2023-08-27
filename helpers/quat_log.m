function r3 = quat_log( quat )

    assert( all( size( quat ) == [ 1,5 ] ) || all( size( quat ) == [ 4,1 ] ) );
    assert( abs( dot( quat, quat ) - 1 ) <= 1e-7 );

    u = quat( 2:4 );
    v = quat( 1 );
    
    if norm( u ) == 0
        r3 = zeros( 3, 1 );
    else
        r3 = acos( v ) * u/norm( u );
    end

end

