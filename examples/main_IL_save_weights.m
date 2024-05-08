% [Title]     Saving the Weights for Imitation Learning
% [Author]    Moses Chong-ook Nah
% [Email]     mosesnah@mit.edu
% [Update]    At 2023.11.15

%% (--) Initialization
clear; close all; clc;
fig_config( 'fontSize', 20, 'markerSize', 10 )

%% =====================================================================
%% (1-) Learning the weights for the Imitation Learning, Discrete Movements
%%  -- (1A) Define Trajectory and its time-trajectory

% Define the time symbol for Analytical Trajectory
syms t_sym

% The three types of trajectory that we aim to learn
% [1] cosine
% [2] radial
% [3] minjerk
name_traj = 'minjerk';

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

    case 'minjerk'

        px = 10 * tn^3 - 15 * tn^4 + 6 * tn^5;
        py = 10 * tn^3 - 15 * tn^4 + 6 * tn^5;
        pz = 10 * tn^3 - 15 * tn^4 + 6 * tn^5;

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
alpha_z = 100.0;
beta_z  = 0.25 * alpha_z;
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
data.weight  = w_arr( 1:2, : );

% Goal Location, where initial position is automatically zero
data.goal =  g_d( 1:2 );
data.z0   = dp_data( 1:2, 1 )/tau;

save( ['learned_parameters/', name_traj ,'.mat' ], 'data' );



%% =====================================================================
%% (2-) Learning the weights for the Imitation Learning, Rhythmic Movements
%%  -- (2A) Define Trajectory and find the best-fit weights

% Define the time symbol for Analytical Trajectory
syms t_sym

% The three types of trajectory that we aim to learn
% [1] circle
% [2] heart
% [3] eight
name_traj = 'eight';

% Period 
Tp = 1;       

switch name_traj

    case 'circle'
        r = 0.5;
        t = 2*pi/Tp * t_sym;
        x = r * cos( t );
        y = r * sin( t );

    case 'heart'
        t = 2*pi/Tp * t_sym;
        x = 16 * sin( t )^3;
        y = 13 * cos( t ) - 5 * cos( 2*t )- 2 * cos( 3*t ) - cos( 4*t );

    case 'eight'
        t = 2*pi/Tp * t_sym;
        x = 2 * sin( t );
        y = 2 * sin( t ) * cos( t );

    otherwise
        error( 'Wrong input: %s', name_traj );  
end

  p_sym = [ x;y ]; 
 dp_sym = diff(  p_sym, t_sym );
ddp_sym = diff( dp_sym, t_sym );

  p_func = matlabFunction(   p_sym );
 dp_func = matlabFunction(  dp_sym );
ddp_func = matlabFunction( ddp_sym );

dt = 1e-4;
t_arr = 0:dt:Tp; 

% Since it is a cyclic input, we do not need the final value
  t_arr = t_arr( 1:end-1 );
  p_arr =   p_func( t_arr );
 dp_arr =  dp_func( t_arr );
ddp_arr = ddp_func( t_arr );

% Get the length of the t_arr
P = length( t_arr );

% Get the initial position
p_arr = p_arr - p_arr( :, 1 );

% Get the average position and take it off
goal = mean( p_arr, 2 );

% Parameters of DMP
N = 50;
alpha_s = 1.0;
alpha_z = 100.0;
beta_z  = 0.25 * alpha_z;
tau = Tp/(2*pi);

% Defining the DMPs
cs        = CanonicalSystem( 'rhythmic', tau, alpha_s );
trans_sys = TransformationSystem( alpha_z, beta_z, cs );
fs        = NonlinearForcingTerm( cs, N );

% The Force Array 
% Saving the initial and goal location of the demonstrated trajectory
f_arr = trans_sys.get_desired( p_arr, dp_arr, ddp_arr, goal );

% The phi matrix 
Phi_mat = zeros( P, N );

% Calculate the Phi matrix for Weight Learning
% Note that compared to Ijspeert 2013, this followed Koutras and Doulgeri (2020)
for i = 1 : P 
    t = t_arr( i );
    Phi_mat( i, : ) = fs.calc_multiple_ith( t, 1:N )/ fs.calc_whole_at_t( t );
end

% Learning the weights with Linear Least-square fitting
w_arr = transpose( ( Phi_mat' * Phi_mat )^-1 * Phi_mat' * f_arr' );

%%  -- (2B) Double Check the Results

T   = Tp;
N   = 3000;
t_arr = linspace( 0, T, N+1 ); 

f = figure( ); a = axes( 'parent', f );
hold( a, 'on' ); axis equal

input_arr = fs.calc_forcing_term( t_arr( 1:end-1 ), w_arr, 0, eye( 2 ) ) + alpha_z*beta_z*goal ;
[ p_roll, ~, ~] = trans_sys.rollout( zeros( 2, 1 ), zeros( 2, 1 ), zeros( 2, 1 ), input_arr, 0, t_arr  ); 

plot( a, p_roll( 1, : ), p_roll( 2, : ), 'linewidth', 4 )
plot( a,  p_arr( 1, : ),  p_arr( 2, : ), 'linewidth', 8, 'linestyle', '--' , 'color', 'k')

%%  -- (2C) Saving the Data for future rollout

data = struct;

% Parameters of DMP and Learned Weighted
data.name    = name_traj;
data.type    = 'rhythmic';
data.tau     = tau;
data.alpha_z = alpha_z;
data.beta_z  =  beta_z;
data.weight  = w_arr;
data.p_init  =  p_arr( :, 1 );
data.dp_init =  dp_func( 0 );
data.goal    = goal;

save( ['learned_parameters/rhythmic/', name_traj ,'.mat' ], 'data' );
