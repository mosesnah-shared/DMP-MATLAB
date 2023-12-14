function [ pos_arr, vel_arr, acc_arr ] = min_jerk_traj( xi, xf, D, t_arr, t0i )
% ===========================================================================
% Descriptions
% ------------
%    Minimum Jerk Trajectory in n dimensional space
% 
% Parameters
% ----------
%   (1) xi: initial position, row or column vector
%  
%   (2) xf: final position, row or column vector
% 
%   (3)  D: duration of the movement
%
%   (4)  t_arr: time array of the trajectory
%
%   (5)  t0i: initial time of the movement
% 
% Returns
% -------
%   (1)   x: position 
%
%   (2)  dx: velocity
%
%   (3) ddx: acceleration
%
% ===========================================================================

% Initial and final positions must be either row or column vectors
assert( ( isrow( xi ) && isrow( xf ) ) || ( iscolumn( xi ) && iscolumn( xf ) ) )

% Must be also the same size
assert( length( xi ) == length( xf ) );

% Get the length of both dimension and time array
n  = length(    xi );
Nt = length( t_arr ); 

% The initial time must be smaller than maximum t_arr
assert( t0i <= max( t_arr ) )

% Generate a minimum-jerk trajectory, pos. vel. and acc.
pos_arr = zeros( n, Nt );
vel_arr = zeros( n, Nt );
acc_arr = zeros( n, Nt );

for i = 1 : Nt
    t = t_arr( i );

    tau = ( t-t0i )/D;

    if t <= t0i 
        pos_arr( :, i ) = xi;

    elseif t >= t0i + D
        pos_arr( :, i ) = xf;

    else
        pos_arr( :, i ) = xi + (xf - xi) * ( 10 * tau^3 -  15 * tau^4 +   6 * tau^5 );
        vel_arr( :, i ) =      (xf - xi) * ( 30 * tau^2 -  60 * tau^3 +  30 * tau^4 )/D^1;
        acc_arr( :, i ) =      (xf - xi) * ( 60 * tau^1 - 180 * tau^2 + 120 * tau^3 )/D^2;
    end
end

end