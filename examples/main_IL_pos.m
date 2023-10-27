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
q0i = [ 0.0; 3.0; 2.0 ];
q0f = [ 2.0; 1.0; 1.0 ];

n   = length( q0i );
D   = 2.0;
t0i = 1.0;

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
trans_sys = TransformationSystem( alpha_z, beta_z, cs );
fs        = NonlinearForcingTerm( cs, N );

%% ---- [1B] Learning Weights via Locally Weighted Regression or Least-Square

% Assume that we have P sample points for the minimum jerk trajectory.
P = 100;

% Equal Sampling along time with duration D
t_P = linspace( 0.0, D, P );

% Desired Trajectories
  y_des_arr = zeros( 3, P );
 dy_des_arr = zeros( 3, P );
ddy_des_arr = zeros( 3, P );

for i = 1 : 3
    for j = 1 : P
        [ y_des, dy_des, ddy_des ] = min_jerk_traj( t_P( j ), q0i, q0f, D, 0 );
          y_des_arr( :, j ) =   y_des;
         dy_des_arr( :, j ) =  dy_des;
        ddy_des_arr( :, j ) = ddy_des;
    end
end

f_arr = trans_sys.get_desired( y_des_arr, dy_des_arr, ddy_des_arr, g ); 

% Learning the weights with least-square method 
w_arr_LSS = zeros( n, N );

for k = 1 : n
    phi_mat = zeros( P, N );
    
    for i = 1 : P
        phi_sum = fs.calc_whole_at_t( t_P( i ) );
        for j = 1 : N
            phi_mat( i, j ) = fs.calc_ith( t_P( i ), j ) / phi_sum * ( g( k ) - y0( k ) ) * cs.calc( t_P( i ) );
        end
    end
    
    % Get w_arr with Least square solution
    w_arr_LSS( k, : ) = transpose( ( phi_mat' * phi_mat )^(-1) * phi_mat' * f_arr( k, : )' );

end

%% ---- [1C] Rollout Generating a Full trajectory with the Transformation System



%% ---- [1D] Saving the Weights and Parameters of the Learned Trajectory 

% All the necessary data 
data = struct;

data.alpha_z = alpha_z;
data.beta_z  = beta_z;

data.cs        = cs;
data.trans_sys = trans_sys;
data.fs        = fs;
data.g         = g;

data.w_arr_LSS = w_arr_LSS;
data.w_arr_LWR = w_arr_LWR;

% save( './learned_parameters/min_jerk_traj' , 'data' );
