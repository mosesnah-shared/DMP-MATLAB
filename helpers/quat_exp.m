function quat_new = quat_exp( w )

    assert( all( size( w ) == [ 1,3 ] ) || all( size( w ) == [ 3,1 ] ) );
    w_norm = norm( w );
    
    if w_norm <=1e-9
       quat_new = [1,0,0,0]' ;
    else
       quat_new = [ cos( w_norm ); sin( w_norm ) * w / w_norm ];
    end

end

