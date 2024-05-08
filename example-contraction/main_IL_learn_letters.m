% [Title]     Saving the Weights for Imitation Learning, Import DAta
% [Author]    Moses Chong-ook Nah
% [Email]     mosesnah@mit.edu
% [Update]    At 2023.11.15

%% (--) Initialization
clear; close all; clc;
fig_config( 'fontSize', 20, 'markerSize', 10 )

%% =====================================================================
%% (1-) Learning the Weights of the Alphabet
%%  -- (1A) Import Data

% Import the data
traj_name = 'hill_shape';
load( [ './alphabets/', traj_name, '.mat' ] );

% Time, position, velocity and acceleration of the data
  t_data = data.t_arr;
  p_data = data.p_data;
 dp_data = data.dp_data;
ddp_data = data.ddp_data;

% Duration and number of samples
tau =    max( t_data );
P   = length( t_data );

% Initialize the position to zero
p_data = p_data - p_data( :, 1 );

%%  -- (1B) Imitation Learning

% Parameters for DMP and Imitation Learning
N       = 50;
alpha_s = 1.0;
alpha_z = 100;
beta_z  = 0.25 * alpha_z;

% Defining the DMPs
cs        = CanonicalSystem( 'discrete', tau, alpha_s );
trans_sys = TransformationSystem( alpha_z, beta_z, cs );
fs        = NonlinearForcingTerm( cs, N );

% Calculating the Force array
g_d   = p_data( :, end );
B_mat = trans_sys.get_desired( p_data, dp_data, ddp_data, g_d );

% The phi matrix 
A_mat = zeros( N, P );

% Calculate the Phi matrix for Weight Learning
% Note that compared to Ijspeert 2013, this followed Koutras and Doulgeri (2020)
for i = 1 : P 
    t = t_data( i );
    A_mat( :, i ) = fs.calc_multiple_ith( t, 1:N )/ fs.calc_whole_at_t( t ) * cs.calc( t );
end

% Learning the weights with Linear Least-square fitting
w_arr = B_mat * A_mat' * inv( A_mat * A_mat' );

% Rollout the trajectory with the weight array and compare its output
t0i   = 0.3;
T     = tau+3;
dt    = 1e-4;
t_arr = 0:dt:T;

% New initial condition
z0_new = dp_data( :, 1 )*tau;
g_new  = p_data( :, end );

f = figure( ); a = axes( 'parent', f );
axis equal; hold( a, 'on' );

input_arr = fs.calc_forcing_term( t_arr( 1:end-1 ), w_arr, t0i, eye( 2 ) );
[ p_arr, z_arr, dp_arr ] = trans_sys.rollout( zeros( 2, 1 ), z0_new, g_new, input_arr, t0i, t_arr  ); 

% Plotting the Line
plot( a,  p_arr( 1, : ),  p_arr( 2, : ), 'linewidth', 5, 'color', 'k' )
plot( a, p_data( 1, : ), p_data( 2, : ), 'linewidth', 5, 'color', 'g' )

% The initial and final position
scatter( a, p_arr( 1,   1 ), p_arr( 2,   1 ), 'filled', 'd', 'markerfacecolor', 'white', 'markeredgecolor', [0.0000, 0.4470, 0.7410], 'linewidth', 5 )
scatter( a, p_arr( 1, end ), p_arr( 2, end ), 'filled', 'o', 'markerfacecolor', 'white', 'markeredgecolor', [0.8500, 0.3250, 0.0980], 'linewidth', 5 )

% The initial and final position
% scatter( a, p_arr( 1,   1 ), p_arr( 2,   1 ), 'filled', 'd', 'markerfacecolor', 'white', 'markeredgecolor', [0.0000, 0.4470, 0.7410], 'linewidth', 5 )
scatter( a, g_new( 1 ), g_new( 2 ), 'filled', 'o', 'markerfacecolor', 'white', 'markeredgecolor', [0.8500, 0.3250, 0.0980], 'linewidth', 5 )


%%  -- (1C) Saving the Data for future rollout

data = struct;

% Parameters of DMP and Learned Weighted
data.name    = traj_name;
data.tau     = tau;
data.alpha_s = alpha_s;
data.alpha_z = alpha_z;
data.beta_z  =  beta_z;
data.weight  = w_arr;

% Goal Location, where initial position is automatically zero
data.goal =  g_d;
data.z0   = dp_data( :, 1 )/tau;

save( ['learned_parameters/discrete/', traj_name ,'_loose.mat' ], 'data' );
