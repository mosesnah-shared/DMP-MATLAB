function quat = axang_to_quat( theta, axis )
    
    quat = zeros( 1, 4 );
    quat(   1  ) = cos( theta/2 );
    quat( 2:4  ) = sin( theta/2 ) * axis;
       
end

