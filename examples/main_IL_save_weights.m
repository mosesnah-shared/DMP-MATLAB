% [Title]     Saving the Weights for Imitation Learning
% [Author]    Moses Chong-ook Nah
% [Email]     mosesnah@mit.edu
% [Update]    At 2023.11.15

%% (--) Initialization
clear; close all; clc;
fig_config( 'fontSize', 20, 'markerSize', 10 )

%% =====================================================================
%% (1-) Learning the weights for the Imitation Learning, Analytical Form
%%  -- (1A) Define Trajectory and its time-trajectory

% Define the time symbol for Analytical Trajectory
syms t_sym

% The three types of trajectory that we aim to learn
% [1] cosine
% [2] radial
name_traj = 'radial';

% Duration and Normalized Time        
D  = 3.0;   
tn = t_sym/D;       

switch name_traj

    case 'cosine'
        % Parameters
        l  = 1.0;

        px = 0.0;
        py =   l * ( 10 * tn^3 - 15 * tn^4 + 6 * tn^5 );
        pz = 0.5 * ( cos( 2*2*pi/l * ( l * tn ) ) - 1 );

    case 'radial'
        % Parameters 
        r0i   = 0.45;
        r0f   = 0.15;
        rt    = r0i + ( r0f - r0i ) * ( -2*tn^3 + 3*tn^2 );
        theta = 2 * pi * ( -tn^4 + 2*tn^2 );

        px = 0;
        py = rt * cos( theta );
        pz = rt * sin( theta );

    otherwise
        error( 'Wrong input: %s', name_traj )
end
    
% Velocity and Acceleration
 dpx = diff(  px, t_sym );  dpy = diff(  py, t_sym );  dpz = diff(  pz, t_sym );
ddpx = diff( dpx, t_sym ); ddpy = diff( dpy, t_sym ); ddpz = diff( dpz, t_sym );

% Saving as MATLAB functions from symbolic form
  p_func = matlabFunction( [   px;   py;   pz ] );
 dp_func = matlabFunction( [  dpx,  dpy,  dpz ] );
ddp_func = matlabFunction( [ ddpx, ddpy, ddpz ] );

% Generating the actual data from the symbolic form 
% These data will be used for Imitation Learning
dt    = 1e-2;
t_arr = 0:dt:D;
P     = length( t_arr );

% Initialization of the data point for position, velocity and acceleration
  p_data = zeros( 3, P );
 dp_data = zeros( 3, P );
ddp_data = zeros( 3, P );

% For the trajectory learning, we simply
for i = 1 : P
    t = t_arr( i );
      p_data( :, i ) =   p_func( t );
     dp_data( :, i ) =  dp_func( t );
    ddp_data( :, i ) = ddp_func( t );
end

% Adding an offset to make the initial position 0
p_data = p_data - p_data( :, 1 );

% Plotting the trajectory
f = figure( ); a = axes( 'parent', f );
plot3( a, p_data( 1, : ), p_data( 2, : ), p_data( 3, : ), 'linewidth', 5, 'color', 'black' )
axis equal 
view( 3 );

%%  -- (1B) Learning the weights of the demonstrated trajectory

% For Imitation Learning, one should define the number of basis function
N  = 50;

% Parameters of DMP
alpha_s = 1.0;
alpha_z = 1000.0;
beta_z  = 0.5 * alpha_z;
tau = D;

% Defining the DMPs
cs        = CanonicalSystem( 'discrete', D, alpha_s );
trans_sys = TransformationSystem( alpha_z, beta_z, cs );
fs        = NonlinearForcingTerm( cs, N );

% The Force Array 
% Saving the initial and goal location of the demonstrated trajectory
pi_d  = p_data( :,   1 );
g_d   = p_data( :, end );
f_arr = trans_sys.get_desired( p_data, dp_data, ddp_data, g_d );

% The phi matrix 
Phi_mat = zeros( P, N );

% Calculate the Phi matrix for Weight Learning
% Note that compared to Ijspeert 2013, this followed Koutras and Doulgeri (2020)
for i = 1 : P 
    t = t_arr( i );
    Phi_mat( i, : ) = fs.calc_multiple_ith( t, 1:N )/ fs.calc_whole_at_t( t ) * cs.calc( t );
end

% Learning the weights with Linear Least-square fitting
w_arr = transpose( ( Phi_mat' * Phi_mat )^-1 * Phi_mat' * f_arr' );

% Rollout with the weight array 
t0i = 1.0;
T   = 5.0;
N   = 3000;
t_arr = linspace( 0, T, N+1 ); 

% New initial condition
z0_new = dp_data( :, 1 )/D;

f = figure( ); a = axes( 'parent', f );
view( 3 ); axis equal
hold( a, 'on' );

rot1 = roty( 30 );
scl  = 1.0;

for angle = 0:120:360
    g_new  = scl * rotz( angle ) * rot1 *  g_d;
    input_arr = fs.calc_forcing_term( t_arr( 1:end-1 ), w_arr, t0i, scl * rotz( angle ) * rot1 );
    [ p_arr, z_arr, dp_arr ] = trans_sys.rollout( zeros( 3, 1 ), z0_new, g_new, input_arr, t0i, t_arr  ); 

    % Plotting the Line
    plot3( a, p_arr( 1, : ), p_arr( 2, : ), p_arr( 3, : ), 'linewidth', 5, 'color', 'k' )

    % The initial and final position
    scatter3( a, p_arr( 1,   1 ), p_arr( 2,   1 ), p_arr( 3,   1 ), 'filled', 'd', 'markerfacecolor', 'white', 'markeredgecolor', [0.0000, 0.4470, 0.7410], 'linewidth', 5 )
    scatter3( a, p_arr( 1, end ), p_arr( 2, end ), p_arr( 3, end ), 'filled', 'o', 'markerfacecolor', 'white', 'markeredgecolor', [0.8500, 0.3250, 0.0980], 'linewidth', 5 )
end

%%  -- (1C) Saving the Data for future rollout

data = struct;

% Parameters of DMP and Learned Weighted
data.name    = name_traj;
data.tau     = tau;
data.alpha_s = alpha_s;
data.alpha_z = alpha_z;
data.beta_z  =  beta_z;
data.weight  = w_arr;

% Goal Location, where initial position is automatically zero
data.goal =  g_d;
data.z0   = dp_data( :, 1 )/tau;

save( ['learned_parameters/', name_traj ,'.mat' ], 'data' );
