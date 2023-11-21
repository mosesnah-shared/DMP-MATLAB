function new_quat = quat_power( quat, t )
    
    % Check unit quaternion
    assert( is_unit_quat( quat ) );

    % Get the axis angle
    [ theta, axis ] = quat_to_axang( quat );

    new_quat = axang_to_quat( theta/2 * t, axis );
       
end

