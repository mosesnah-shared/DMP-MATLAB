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

% Generating Minimum-jerk trajectory
qid = [ 0.0; 0.0 ];         % Initial Posture, Demonstrated
qgd = [ 1.0; 1.0 ];         %   Final Posture, Demonstrated 
D   = 1.0;                  % Duration
n   = length( qid );        % Number of Degrees of Freedom

% Parameters for the Canonical System, Discrete
tau     =    D;
alpha_s =  1.0;

% The number of Basis Function for the Nonlinear Forcing Term
N = 25;

% Parameters of the Transformation System
alpha_z   = 50.0;
beta_z    = 0.25 * alpha_z;

% DMP, three primitives
cs        = CanonicalSystem( 'discrete', tau, alpha_s );
trans_sys = TransformationSystem( alpha_z, beta_z, cs );
fs        = NonlinearForcingTerm( cs, N );

% P sample points for the minimum jerk trajectory.
P = 100;

% Equal Sampling along time with duration D
t_P = linspace( 0.0, D, P );

[ q_des, dq_des, ddq_des ] = min_jerk_traj( qid, qgd, D, t_P, 0 );

% Calculating the required Nonlinear Forcing Term
B = trans_sys.get_desired( q_des, dq_des, ddq_des, qgd );

% The A matrix 
A = zeros( N, P );

% Interating along the sample points
for i = 1 : P 
    t = t_P( i );
    A( :, i ) = fs.calc_multiple_ith( t, 1:N )/ fs.calc_whole_at_t( t ) * cs.calc( t );
end

% Scaling matrix, classical DMP
scl = diag( qgd - qid );

% Linear Least-square
W = B * A' * ( A * A' )^-1;

% Rollout with the weight array 
t0i   = 0.0;
T     = 1.0;
Nt    = 3000;
t_arr_dis = linspace( 0, T, Nt+1 );

[ q_des_plot , ~, ~] = min_jerk_traj( qid, qgd, D, t_arr_dis, t0i );

input_arr = fs.calc_forcing_term( t_arr_dis( 1:end-1 ), W, t0i, diag( qgd - qid ) );

[ y_arr_dis, ~, ~ ] = trans_sys.rollout( q_des( :, 1 ), dq_des( :, 1 ), q_des( :, end ), input_arr, t0i, t_arr_dis  );

f = figure( );
a = subplot( 2, 1, 1 );
hold on;
% First Joint
plot( t_arr_dis, y_arr_dis( 1, : ), 'linewidth', 5, 'color', [0, 0.4470, 0.7410], 'linestyle', '-' );
plot( t_arr_dis, q_des_plot( 1, : ), 'linewidth', 8, 'color', 'k', 'linestyle', '--' );
set( a, 'fontsize', 25, 'xtick', 0.0:0.5:1.5, 'xticklabel', {},'ylim', [ qid(1)-0.2, qgd(1)+0.2 ] );
ylabel( a, 'Joint 1 (rad)', 'fontsize', 35 )

a = subplot( 2, 1, 2 );
% Second Joint
hold on;
plot( t_arr_dis, y_arr_dis( 2, : ), 'linewidth', 5, 'color', [0, 0.4470, 0.7410], 'linestyle', '-' );
plot( t_arr_dis, q_des_plot( 2, : ), 'linewidth', 8, 'color', 'k', 'linestyle', '--' );
set( a, 'fontsize', 25, 'xtick', 0.0:0.5:1.5, 'ylim', [ qid(2)-0.2, qgd(2)+0.2 ] );
xlabel( a, 'Time (sec)', 'fontsize', 35 )
ylabel( a, 'Joint 2 (rad)', 'fontsize', 35 )

%%  -- (1B) For Rhythmic Movement

% Generating rhythmic movement, cosine and sine
r0  = 0.5;
Tp  = 2.0;
tau = Tp/(2*pi);
w   = 1/tau;        % Angular velocity of cosine/sin

% The number of Basis Function for the Nonlinear Forcing Term
% The algorithm is sensitive to the choice of N!! Care is required.
N = 25;

% Parameters of the Transformation System
alpha_z   = 50.0;
beta_z    = 0.25 * alpha_z;

P   = 100;
t_P = linspace( 0, Tp, P );
  q_des =        r0 * [  cos( w * t_P );  sin( w * t_P ) ];
 dq_des =  w^1 * r0 * [ -sin( w * t_P );  cos( w * t_P ) ];
ddq_des =  w^2 * r0 * [ -cos( w * t_P ); -sin( w * t_P ) ];

% DMP
cs        = CanonicalSystem( 'rhythmic', tau, 1.0 );
trans_sys = TransformationSystem( alpha_z, beta_z, cs );
fs        = NonlinearForcingTerm( cs, N );

% Calculating the required Nonlinear Forcing Term
% Goal location is zero 
B = trans_sys.get_desired( q_des, dq_des, ddq_des, zeros( 2, 1 ) );

% The A matrix 
A = zeros( N, P );

% Interating along the sample points
for i = 1:P 
    t = t_P( i );
    A( :, i ) = fs.calc_multiple_ith( t, 1:N )/ fs.calc_whole_at_t( t );
end

% Linear Least-square
W = B * A' * ( A * A' )^-1;

T  = Tp*2;
Nt = 5000;
t_arr_rhy = linspace( 0, T, Nt+1 );

input_arr = fs.calc_forcing_term( t_arr_rhy( 1:end-1 ), W, 0, eye( 2 ) );
[ y_arr_rhy, ~, dy_arr_dis ] = trans_sys.rollout( q_des( :, 1 ), tau*dq_des( :, 1 ), zeros( 2, 1 ), input_arr, 0, t_arr_rhy  );

f = figure( );
a = subplot( 2, 1, 1 );
% First Joint
hold on;
plot( t_arr_rhy, y_arr_rhy( 1, : ), 'linewidth', 5, 'color', [0, 0.4470, 0.7410], 'linestyle', '-' );
plot( t_arr_rhy, r0 * cos( w*t_arr_rhy ), 'linewidth', 8, 'color', 'k', 'linestyle', '--' );
set( a, 'fontsize', 25, 'ylim', [-r0-0.1, r0+0.1]  );
ylabel( a, 'Joint 1 (rad)', 'fontsize', 35 )

a = subplot( 2, 1, 2 );
% Second Joint
hold on;
plot( t_arr_rhy, y_arr_rhy( 2, : ), 'linewidth', 5, 'color', [0, 0.4470, 0.7410], 'linestyle', '-' );
plot( t_arr_rhy, r0 * sin( w*t_arr_rhy ), 'linewidth', 8, 'color', 'k', 'linestyle', '--' );
set( a, 'fontsize', 25, 'ylim', [-r0-0.1, r0+0.1] );
xlabel( a, 'Time (sec)', 'fontsize', 35 )
ylabel( a, 'Joint 2 (rad)', 'fontsize', 35 )

%%  -- (1C) For both Discrete and Rhythmic Movement

f = figure( ); 

a = subplot( 2, 2, 1 );
hold on;
plot( a, t_arr_dis,  y_arr_dis( 1, : ), 'linewidth', 5, 'color',[0, 0.4470, 0.7410], 'linestyle', '-' );
plot( a, t_arr_dis, q_des_plot( 1, : ), 'linewidth', 8, 'color', 'k', 'linestyle', '--' );
set( a, 'ylim', [ -0.1, 1.1], 'ytick', 0:0.5:1.0, 'fontsize', 25, 'xtick', 0:0.5:1.0, 'xticklabel', {} )
ylabel( a, 'Joint 1 (rad)', 'fontsize', 30 )
title( a, 'Discrete', 'fontsize', 35 )

a = subplot( 2, 2, 3 );
hold on;
plot( a, t_arr_dis,  y_arr_dis( 2, : ), 'linewidth', 5, 'color',[0, 0.4470, 0.7410], 'linestyle', '-' );
plot( a, t_arr_dis, q_des_plot( 2, : ), 'linewidth', 8, 'color', 'k', 'linestyle', '--' );
set( a, 'ylim', [ -0.1, 1.1], 'ytick', 0:0.5:1.0, 'fontsize', 25, 'xtick', 0:0.5:1.0 )
ylabel( 'Joint 2 (rad)', 'fontsize', 30 )
xlabel( 'Time (sec)', 'fontsize', 30 )

a = subplot( 2, 2, 2 );
hold on;
plot( a, t_arr_rhy,       y_arr_rhy( 1, : ), 'linewidth', 5, 'color',[0, 0.4470, 0.7410], 'linestyle', '-' );
plot( a, t_arr_rhy, r0 * cos( w*t_arr_rhy ), 'linewidth', 8, 'color', 'k', 'linestyle', '--' );
set( a, 'ylim', [ -r0-0.1, r0+0.1], 'xtick', 0:4, 'fontsize', 25, 'ytick', -0.5:0.5:0.5, 'xticklabel', {} )
title( a, 'Rhythmic', 'fontsize', 35 )

a = subplot( 2, 2, 4 );
hold on;
plot( a, t_arr_rhy,       y_arr_rhy( 2, : ), 'linewidth', 5, 'color',[0, 0.4470, 0.7410], 'linestyle', '-' );
plot( a, t_arr_rhy, r0 * sin( w*t_arr_rhy ), 'linewidth', 8, 'color', 'k', 'linestyle', '--' );
set( a, 'ylim', [ -r0-0.1, r0+0.1], 'xtick', 0:4, 'fontsize', 25, 'ytick', -0.5:0.5:0.5 )
xlabel( 'Time (sec)', 'fontsize', 30 )

fig_save( f, 'ThesisImages/images/joint_space_imit_learn' )

%% =======================================================
%% (2-) Imitation Learning for Task-space, Position, Section 4.4.2
%%  -- (2A) For Discrete Movement

% Load the M data
data = load( 'learned_parameters/M.mat' ); data = data.data;

% System parameters
tau = 1.0;
alpha_s = data.alpha_s;
alpha_z = data.alpha_z;
 beta_z =  data.beta_z;

W = data.weight;
N = size( W, 2 );

% DMP, three primitives
cs        = CanonicalSystem( 'discrete', tau, alpha_s );
trans_sys = TransformationSystem( alpha_z, beta_z, cs );
fs        = NonlinearForcingTerm( cs, N );

pg = data.goal;
pi = zeros( 2, 1 );

t0i = 0.0;

% Nt is the sample point
Nt = 5000;
t_arr_dis = linspace( 0, tau, Nt+1 );

input_arr = fs.calc_forcing_term( t_arr_dis( 1:end-1 ), W, t0i, eye( 2 ) );
[ y_arr_dis, ~, ~ ] = trans_sys.rollout( pi, data.z0/tau, pg, input_arr, t0i, t_arr_dis  );

plot( y_arr_dis( 1, : ),y_arr_dis( 2, : ))
%%  -- (2B) For Rhythmic Movement


%%  -- (2C) For Discrete and Rhythmic Movement



% For discrete and rhythmic movement