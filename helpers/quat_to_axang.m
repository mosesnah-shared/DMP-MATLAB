function [ theta, axis ] = quat_to_axang( quat )
    
    % Check unit quaternion
    assert( is_unit_quat( quat ) );

    qw = quat_real( quat );
    qv = quat_imag( quat );

    if norm( qv ) <= 1e-9
        axis = zeros( 1, 3 );
    else
        axis = qv/norm( qv );
    end

    % Get theta
    theta = 2*atan2( norm( qv ), qw );
       
end

