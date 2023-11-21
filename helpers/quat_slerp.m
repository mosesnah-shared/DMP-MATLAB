function quat_new = quat_slerp( quat1, quat2, t_arr, t0i, t0f )
% Interpolation between 
assert( iscolumn( t_arr ) || isrow( t_arr ) );
assert( max( t_arr ) >= t0f && t0i >= 0 && t0i <= t0f );

% Zero offset
t_arr = t_arr - t_arr( 1 );

Nt = length( t_arr );
quat_new = zeros( 4, Nt );

for i = 1 : Nt
    t = t_arr( i );

    if t <= t0i 
        quat_new( :, i ) = quat1;
    elseif t >= t0f
        quat_new( :, i ) = quat2;
    else
        tau = ( t-t0i )/( t0f - t0i );
        
        % Dot product
        tmp = acos( dot( quat1, quat2 ) );
        
        quat_new( :, i ) = quat1 * sin( (1-tau)*tmp )/sin( tmp ) + quat2 * sin( tau*tmp ) /sin(tmp);
    end
end

