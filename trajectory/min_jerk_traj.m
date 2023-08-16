function [ x, dx, ddx ] = min_jerk_traj( t, x0i, x0f, D, ti )
    
    if t <= ti 
         x = x0i; 
        dx =   0;
       ddx =   0;
       
    elseif ti < t && t <= ti + D
         tau = ( t - ti )/D;
           x =     x0i + ( x0f - x0i ) * ( 10 * tau^3 -  15 * tau^4 +   6 * tau^5 );
          dx = 1.0/D   * ( x0f - x0i ) * ( 30 * tau^2 -  60 * tau^3 +  30 * tau^4 );
         ddx = 1.0/D^2 * ( x0f - x0i ) * ( 60 * tau   - 180 * tau^2 + 120 * tau^3 );
    else
         x = x0f; 
        dx =   0;
       ddx =   0;
    end

end

