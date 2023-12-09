% [Title]     Generating Images for Chapter 4
% [Author]    Moses Chong-ook Nah
% [Email]     mosesnah@mit.edu
% [Update]    At 2023.12.05

%% (--) Initialization
clear; close all; clc;
fig_config( 'fontSize', 20, 'markerSize', 10 )

%% =======================================================
%% (1-) Imitation Learning for Joint-space, Section 4.4.1
%%  -- (1A) For Discrete Movement

% The trajectory we aim to imitate is the Minimum jerk Trajectory
% The parameters of the trajectory include:
%   (1) Initial position (y0d) (column vec.)
%   (2)   Final position (gd) (column vec.)
%   (3) Duration (D) 
%   (4) Starting time (t0i)
% For (1), (2) the g subscript is to emphasize the initial/final
% position of a "d"emonstrated trajectory

qid = [ 0.0; 0.0 ];
qgd = [ 1.0; 1.0 ];
D   = 1.0;
n   = length( qid );

% Parameters for the Canonical System
tau     =    D;
alpha_s =  1.0;

% The number of Basis Function for the Nonlinear Forcing Term
N = 50;

% Parameters of the Transformation System
alpha_z   = 50.0;
beta_z    = 0.25 * alpha_z;

cs        = CanonicalSystem( 'discrete', tau, alpha_s );
trans_sys = TransformationSystem( alpha_z, beta_z, cs );
fs        = NonlinearForcingTerm( cs, N );

% P sample points for the minimum jerk trajectory.
P = 100;

% Equal Sampling along time with duration D
t_P = linspace( 0.0, D, P );

% Desired Trajectories, position, velocity and acceleration.
  q_des = zeros( n, P );
 dq_des = zeros( n, P );
ddq_des = zeros( n, P );

for i = 1 : P
    [ q, dq, ddq ] = min_jerk_traj( t_P( i ), qid, qgd, D, 0 );
      q_des( :, i ) =   q;
     dq_des( :, i ) =  dq;
    ddq_des( :, i ) = ddq;
end

% Calculating the required Nonlinear Forcing Term
B = trans_sys.get_desired( q_des, dq_des, ddq_des, qgd );

% The A matrix 
A = zeros( N, P );

% Interating along the sample points
for i = 1 : P 
    t = t_P( i );
    A( :, i ) = fs.calc_multiple_ith( t, 1:N )/ fs.calc_whole_at_t( t ) * cs.calc( t );
end

% Scaling matrix
% There are two options, from Koutras and Doulgeri (2020) and the Ijspeert et al. (2013).
% For this time, we use Ijspeert et al. (2013).
scl = diag( qgd - qid );

W = B * A' * ( A * A' )^-1;

% Rollout with the weight array 
t0i   = 0.0;
T     = 1.0;
Nt    = 3000;
t_arr = linspace( 0, 1, Nt+1 ); 

input_arr = fs.calc_forcing_term( t_arr( 1:end-1 ), W, t0i, diag( qgd - qid ) );

[ y_arr, z_arr, dy_arr ] = trans_sys.rollout( q_des( :, 1 ), dq_des( :, 1 ), q_des( :, end ), input_arr, t0i, t_arr  );

f = figure( ); a = axes( 'parent', f );
hold( a, 'on' )
plot( a, t_arr,  y_arr', 'linewidth', 3 )
plot( a, t_arr, dy_arr', 'linewidth', 3 )

%%  -- (1A) For Rhythmic Movement
