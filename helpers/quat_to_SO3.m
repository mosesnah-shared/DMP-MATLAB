function R = quat_to_SO3( quat )
    % First check whether it is unit quaternion 
    assert( all( size( quat ) == [ 1, 4 ] ) || all( size( quat ) == [ 4, 1 ] ) );
    assert( abs( dot( quat, quat ) - 1 ) <= 1e-7 );

    % get the 1st and 2~4th element
    qw   = quat( 1   );
    qvec = quat( 2:4 );

    qvec_mat = [         0, -qvec( 3 ),  qvec( 2 );
                 qvec( 3 ),          0, -qvec( 1 );
                -qvec( 2 ),  qvec( 1 ),          0];

    R = eye( 3 ) + 2 * qw * qvec_mat + 2 * qvec_mat* qvec_mat;

end

