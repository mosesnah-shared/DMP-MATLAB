function is_check = is_skewsym( mat )
    assert( all( size( mat ) == [ 3, 3 ] ) );

    err = norm( mat + mat' );

    if err <= 1e-5
        is_check = true;
    else
        is_check = false;
    end
    
end