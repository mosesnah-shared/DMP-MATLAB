%% ==================================================================
%% [Title] Imitation Learning, in 3D Space
% Author: Moses Chong-ook Nah
%  Email: mosesnah@mit.edu
%   Date: 2023.08.15
%% ==================================================================

%% [0A] Initialization
clear; close all; clc;
fig_config( 'fontSize', 20, 'markerSize', 10 )

%% [1-] Imitation Learning of Minimum-jerk Trajectory
%% ---- [1A] Parameter Initialization

% The trajectory we aim to imitate is the Minimum jerk Trajectory
% The initial (q0i), final posture (q0f), duration (D), starting time (t0i) of the trajectory
q0i = 0.0;
q0f = 1.0;
D   = 1.0;
t0i = 0.0;

% For Imitation Learning, one should define the number of basis function
N  = 50;

% Parameters of the 3 DMPs
alpha_z = 10.0;
alpha_s = 1.0;
beta_z  = 1/4 * alpha_z;
tau     = D;
g       = q0f;
y0      = q0i;
z0      = 0;

cs        = CanonicalSystem( 'discrete', tau, alpha_s );
trans_sys = TransformationSystem( alpha_z, beta_z, tau, y0, z0 );
fs        = NonlinearForcingTerm( cs, N );

%% ---- [1B] Learning Weights via Locally Weighted Regression or Least-Square

% Assume that we have P sample points for the minimum jerk trajectory.
P = 100;

% Equal Sampling along time with duration D
t_P = linspace( 0.0, D, P );

% Desired Trajectories
  y_des_arr = zeros( 1, P );
 dy_des_arr = zeros( 1, P );
ddy_des_arr = zeros( 1, P );

for j = 1 : P
    [ y_des, dy_des, ddy_des ] = min_jerk_traj( t_P( j ), q0i, q0f, D, 0 );
      y_des_arr( j ) =   y_des;
     dy_des_arr( j ) =  dy_des;
    ddy_des_arr( j ) = ddy_des;
end

% Learning the N weights based on Locally Weighted Regression
w_arr_LWR = zeros( 1, N );

% Calculating a_arr and f_arr that can be calculated
a_arr = ( g - y0 ) * cs.calc( t_P );
f_arr = trans_sys.get_desired( y_des_arr, dy_des_arr, ddy_des_arr, g );

% Iterating over the arrays
for i = 1: N
    phi_arr = fs.calc_ith( t_P, i );    
    w_arr_LWR( i ) = sum( a_arr .* phi_arr .* f_arr ) / sum( a_arr .* phi_arr .* a_arr );
end

a1 = subplot( 2, 1, 1 );
hold( a1, 'on' );

% Plotting the Desired Trajectory
plot( a1, t_P,   y_des_arr );
plot( a1, t_P,  dy_des_arr );
plot( a1, t_P, ddy_des_arr );
 
a2 = subplot( 2, 1, 2 );
hold( a2, 'on' );
for i = 1 : N
    plot( a2, t_P, fs.calc_whole_weighted_at_t( t_P, w_arr_LWR ) )
end

% One can also learn the weights by simple least-square solution
% It is simply f_arr = A w problem, where
% f_arr is simply the array above. 
phi_mat = zeros( P, N );

for i = 1 : P
    phi_sum = fs.calc_whole_at_t( t_P( i ) );
    for j = 1 : N
        phi_mat( i, j ) = fs.calc_ith( t_P( i ), j ) / phi_sum * ( g - y0 ) * cs.calc( t_P( i ) );
    end
end

% Get w_arr with Least square solution
w_arr_LSS = ( phi_mat' * phi_mat )^(-1) * phi_mat' * f_arr';
w_arr_LSS = w_arr_LSS';

%% ---- [1C] Rollout Generating a Full trajectory with the Transformation System

trans_sys.reset( )

% The time step of the simulation and its number of iteration
dt = 1e-3;
Nt = 2000;

% The total time and its time array
T     = dt * Nt;
t_arr = dt * (0:(Nt-1));

% For plotting and saving the data
y_arr  = zeros( 1, Nt );
z_arr  = zeros( 1, Nt );
dz_arr = zeros( 1, Nt );

y_arr( 1 ) = y0;
z_arr( 1 ) = z0;

t = 0;
f_input_arr = zeros( 1, Nt );

% One can also check spatial temporance from the weight
% for tau = 0.5:0.1:2.0
% for   g = 0.5:0.1:2.0

for i = 0 : (Nt-1)
    
    % Before conducting the movement
    % Maintaining that posture
    if t <= t0i
        y_arr( i + 1 ) = y0;
        z_arr( i + 1 ) = z0;
        
    % During the movement
    elseif t0i <= t
        
        if t<= t0i + D

            % taking off the initial time offset
            t_tmp = t - t0i;

            % Calculating the input from the weights
            % First, check if whole activation value is 0
            phi_sum = fs.calc_whole_at_t( t_tmp );
            
            f_input = 0;

            if phi_sum ~= 0
                f_input = fs.calc_whole_weighted_at_t( t_tmp, w_arr_LWR )/phi_sum;
                f_input = f_input*(g-y0)*cs.calc( t_tmp );
            end
            
        else
            f_input = 0; 
        end

        f_input_arr( i ) = f_input; 
        [ y, z, dy, dz ] = trans_sys.step( g, f_input, dt );
        y_arr(  i + 1 ) = y;
        z_arr(  i + 1 ) = z;
        dz_arr( i + 1 ) = dz;        
        dz_arr( i + 1 ) = dz;
    end

    t = t + dt;
end

subplot( 3, 1, 1 )
hold on
plot( t_arr, y_arr, 'linewidth', 3 )
plot( t_P, y_des_arr, 'linestyle', '--', 'linewidth', 5, 'color', 'k' )
set( gca, 'xlim', [0, T], 'fontsize', 30, 'xticklabel', {} )
ylabel( 'Pos. (-)' )

subplot( 3, 1, 2 )
hold on
plot( t_arr, z_arr, 'linewidth', 3 )
plot(   t_P, dy_des_arr, 'linestyle', '--', 'linewidth', 5, 'color', 'k' )
set( gca, 'xlim', [0, T], 'fontsize', 30, 'xticklabel', {} )
ylabel( 'Vel. (-)' )

subplot( 3, 1, 3 )
hold on
plot( t_arr, dz_arr, 'linewidth', 3 )
plot(   t_P, ddy_des_arr, 'linestyle', '--', 'linewidth', 5, 'color', 'k' )
set( gca, 'xlim', [0, T], 'fontsize', 30, 'xtick', 0:0.5:T)
xlabel( 'Time (sec)' )
ylabel( 'Acc. (-)' )

