% [Title]     Saving the Weights for Imitation Learning, Import DAta
% [Author]    Moses Chong-ook Nah
% [Email]     mosesnah@mit.edu
% [Update]    At 2023.11.15

%% (--) Initialization
clear; close all; clc;
fig_config( 'fontSize', 20, 'markerSize', 10 )

%% =====================================================================
%% (1-) Learning the weights for the Imitation Learning
%%  -- (1A) Import Data

% Import the data
traj_name = 'M';
load( [ './alphabets/', traj_name, '.mat' ] );
    
   t_arr = data.t_arr;
  p_data = data.p_data;
 dp_data = data.dp_data;
ddp_data = data.ddp_data;

D = max( t_arr );
P = size( p_data, 2 );

p_data = p_data - p_data( :, 1 );

%%  -- (1B) Learning the weights of the demonstrated trajectory

% For Imitation Learning, one should define the number of basis function
N  = 50;

% Parameters of DMP
alpha_s = 1.0;
alpha_z = 1000.0;
beta_z  = 0.5 * alpha_z;
tau = 6;

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
t0i = 0.3;
T   = D+t0i;
N   = 0.8*1e+5;
dt  = 1e-4;
t_arr = dt*(0:N); 

% New initial condition
z0_new = dp_data( :, 1 )/D;
g_new = p_data( :, end );

f = figure( ); a = axes( 'parent', f );
view( 3 ); axis equal
hold( a, 'on' );

scl = 0.1;
input_arr = fs.calc_forcing_term( t_arr( 1:end-1 ), w_arr, t0i, scl * eye( 2 ), 'trimmed' );
[ p_arr, z_arr, dp_arr ] = trans_sys.rollout( zeros( 2, 1 ), z0_new, scl*g_new, input_arr, t0i, t_arr  ); 


% Plotting the Line
plot( a, p_arr( 1, : ), p_arr( 2, : ), 'linewidth', 5, 'color', 'k' )

% The initial and final position
scatter( a, p_arr( 1,   1 ), p_arr( 2,   1 ), 'filled', 'd', 'markerfacecolor', 'white', 'markeredgecolor', [0.0000, 0.4470, 0.7410], 'linewidth', 5 )
scatter( a, p_arr( 1, end ), p_arr( 2, end ), 'filled', 'o', 'markerfacecolor', 'white', 'markeredgecolor', [0.8500, 0.3250, 0.0980], 'linewidth', 5 )

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

save( ['learned_parameters/', traj_name ,'.mat' ], 'data' );
