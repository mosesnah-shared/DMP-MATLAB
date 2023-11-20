function val = clip_func( t, t0i, t0f )
% Clip the value

assert( t0f > t0i )

if t <= t0i 
    val = 0;
elseif t >= t0f
    val = 1;
else
    val = (t - t0i)/(t0f-t0i);
end

