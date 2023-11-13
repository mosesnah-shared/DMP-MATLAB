function is_check = is_SO3( R )
    assert( all( size( R ) == [ 3, 3 ] ) );
    assert( norm( R * R' - eye( 3 ) ) <= 1e-9 && det( R ) );
    is_check = true;
end