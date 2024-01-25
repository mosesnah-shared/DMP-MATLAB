% [Title]     Generating Images for Chapter 4
% [Author]    Moses Chong-ook Nah
% [Email]     mosesnah@mit.edu
% [Update]    At 2023.12.05

%% (--) Initialization
clear; close all; clc;
fig_config( 'fontSize', 20, 'markerSize', 10 )
c_blue = [0, 0.4470, 0.7410];

%% =======================================================
%% (1-) Imitation Learning for Joint-space, Section 4.4.1
%%  -- (1A) For Discrete Movement

% For demonstration, Generating Minimum-jerk trajectory
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
scl_mat = diag( qgd - qid );

% Linear Least-square
W = inv( scl_mat ) * B * A' * ( A * A' )^-1;

% Rollout with the weight array 
t0i   = 0.0;
T     = 1.0;
Nt    = 3000;
t_arr_dis = linspace( 0, T, Nt+1 );

[ q_des_plot , ~, ~] = min_jerk_traj( qid, qgd, D, t_arr_dis, t0i );

% The new initial condition
qi = qid;
qg = qgd;

input_arr = fs.calc_forcing_term( t_arr_dis( 1:end-1 ), W, t0i, diag( qg - qi ) );
[ y_arr_dis, ~, ~ ] = trans_sys.rollout( qi, zeros( n, 1 ), qg, input_arr, t0i, t_arr_dis  );

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

%%  -- (1Aa) For Discrete Movement, temporal and scaling variance

f = figure( ); 
a1 = subplot( 2, 2, 1 );
hold on 

a2 = subplot( 2, 2, 3 );
set( a1, 'fontsize', 30, 'xtick', 0:2, 'ylim', [0., 1.5] )
set( a2, 'fontsize', 30, 'xtick', 0:2, 'ylim', [0., 1.5] )
xlabel( a2, '$t$ (s)', 'fontsize', 30 )
hold on 

tau_arr = [ 1.0, 0.5, 1.5 ] * D;

N = length( tau_arr );

% Rollout with the weight array 
t0i   = 0.0;
T     = 2*1.0;
Nt    = 3000;
t_arr_dis = linspace( 0, T, Nt+1 );

lw_arr = [ 8, 5, 5 ];
lsty_arr = {'-', '-', '-'};
c_arr = [ c_blue; 0.8500 0.3250 0.0980; 0.9290 0.6940 0.1250];

for i = 1 : N
    tau = tau_arr( i ) ;
   
    % Same goal and initial condition
    qi = qid;
    qg = qgd;

    lw   = lw_arr( i );
    lsty = lsty_arr{ i };

    cs.tau = tau;
    input_arr = fs.calc_forcing_term( t_arr_dis( 1:end-1 ), W, t0i, diag( qg - qi ) );
    [ y_arr_dis, ~, ~ ] = trans_sys.rollout( qi, zeros( n, 1 ), qg, input_arr, t0i, t_arr_dis  );
    
    plot( a1, t_arr_dis, y_arr_dis( 1, : ), 'linewidth', lw, 'linestyle', lsty, 'color', c_arr( i, : ), 'DisplayName', ['$\tau=$', num2str( tau, '%.1f' ) ] )
    plot( a2, t_arr_dis, y_arr_dis( 2, : ), 'linewidth', lw, 'linestyle', lsty, 'color', c_arr( i, : ) )
end
title( a1, 'Temporal Scaling', 'fontsize', 45 );
cs.tau = D;
legend( a1, 'location', 'southeast' )
a3 = subplot( 2, 2, 2 );
title( a3, 'Spatial Scaling', 'fontsize', 45 );
hold on

a4 = subplot( 2, 2, 4 );
hold on 
xlabel( a4, '$t$ (s)', 'fontsize', 30 )

qg_arr = [ 1.0, 1.0; 0.8, 1.2; 0.9, 1.4 ];

for i = 1 : N
    qg  = qg_arr( i, : )';
    % Same goal and initial condition
    qi = qid;

    input_arr = fs.calc_forcing_term( t_arr_dis( 1:end-1 ), W, t0i, diag( qg - qi ) );
    [ y_arr_dis, ~, ~ ] = trans_sys.rollout( qi, zeros( n, 1 ), qg, input_arr, t0i, t_arr_dis  );
    
    lw   = lw_arr( i );
    lsty = lsty_arr{ i };
    
    plot( a3, t_arr_dis, y_arr_dis( 1, : ), 'linewidth', lw, 'linestyle', lsty, 'color', c_arr( i, : ), 'DisplayName', [ '$q_{g,1}=', num2str( qg( 1 ), '%.1f' ), '$' ] )
    plot( a4, t_arr_dis, y_arr_dis( 2, : ), 'linewidth', lw, 'linestyle', lsty, 'color', c_arr( i, : ), 'DisplayName', [ '$q_{g,2}=', num2str( qg( 2 ), '%.1f' ), '$' ] )
end

legend( a3, 'location', 'southeast' )
legend( a4, 'location', 'southeast' )


set( a1, 'xlim', [ 0, 2 ], 'ylim', [0.0, 1.5], 'fontsize', 30, 'xtick', 0:2 )
set( a2, 'xlim', [ 0, 2 ], 'ylim', [0.0, 1.5], 'fontsize', 30, 'xtick', 0:2 )
set( a3, 'xlim', [ 0, 2 ], 'ylim', [0.0, 1.5], 'fontsize', 30, 'xtick', 0:2 )
set( a4, 'xlim', [ 0, 2 ], 'ylim', [0.0, 1.5], 'fontsize', 30, 'xtick', 0:2 )
ylabel( a1, 'Joint 1 (rad)', 'fontsize', 30 )
ylabel( a2, 'Joint 2 (rad)', 'fontsize', 30 )

fig_save( f, 'ThesisImages/images/joint_space_imit_learn_temp_scl' )


%%  -- (1B) For Rhythmic Movement

% Generating rhythmic movement, cosine and sine
r0  = 0.5;
Tp  = 2.0;
tau = Tp/(2*pi);
w   = 1/tau;        % Angular velocity of cosine/sin
r   = 1.0;

% The number of Basis Function for the Nonlinear Forcing Term
% The algorithm is sensitive to the choice of N!! Care is required.
N = 25;

% Parameters of the Transformation System
alpha_z = 50.0;
beta_z  = 0.25 * alpha_z;

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
B = trans_sys.get_desired( q_des, dq_des, ddq_des, mean( q_des, 2 ) );

% The A matrix 
A = zeros( N, P );

% Interating along the sample points
for i = 1:P 
    t = t_P( i );
    A( :, i ) = r * fs.calc_multiple_ith( t, 1:N )/ fs.calc_whole_at_t( t );
end

% Linear Least-square
W = B * A' * ( A * A' )^-1;

T  = Tp*2;
Nt = 5000;
t_arr_rhy = linspace( 0, T, Nt+1 );

input_arr = fs.calc_forcing_term( t_arr_rhy( 1:end-1 ), W, 0, eye( 2 ) );
[ y_arr_rhy, ~, ~ ] = trans_sys.rollout( q_des( :, 1 ), zeros( 2, 1), zeros( 2, 1 ), input_arr, 0, t_arr_rhy  );

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

%%  -- (1Ba) For Rhythmic Movement, temporal and scaling variance

% Create a subplot for 2 column
f = figure( ); 
a1 = subplot( 2, 2, 1 ); hold on
a2 = subplot( 2, 2, 2 ); hold on
a3 = subplot( 2, 2, 3 ); hold on
a4 = subplot( 2, 2, 4 ); hold on

% The first/second columns are Temporal/Spatial Scaling
title( a1, 'Temporal Scaling', 'fontsize', 45 );
title( a2, 'Spatial Scaling', 'fontsize', 45 );

% For the whole rollout
Tp = 2.0;
T  = Tp*2;
Nt = 5000;
t_arr_rhy = linspace( 0, T, Nt+1 );

% Iterate over the temporal scaling
Tp_arr = [ 2.0, 4.0, 1.0 ];
N = length( Tp_arr );

tau_arr = Tp_arr/(2*pi);


lw_arr = [ 8, 5, 5 ];
lsty_arr = {'-', '-', '-'};
c_arr = [ c_blue; 0.8500 0.3250 0.0980; 0.9290 0.6940 0.1250];

for i = 1 : N
    cs.tau = tau_arr( i );

    lw   = lw_arr( i );
    lsty = lsty_arr{ i };

    input_arr = fs.calc_forcing_term( t_arr_rhy( 1:end-1 ), W, 0, eye( 2 ) );
    [ y_arr_rhy, ~, ~ ] = trans_sys.rollout( q_des( :, 1 ), zeros( 2, 1), zeros( 2, 1 ), input_arr, 0, t_arr_rhy  );
    plot( a1, t_arr_rhy, y_arr_rhy( 1, : ), 'linewidth', lw, 'color', c_arr( i, : ), 'linestyle', lsty );
    plot( a3, t_arr_rhy, y_arr_rhy( 2, : ), 'linewidth', lw, 'color', c_arr( i, : ), 'linestyle', lsty );
end

% We try multiple scaling too 
scl_arr = [ 1.0 , 0.5, 1.5 ];
for i = 1 : N
    lw   = lw_arr( i );
    lsty = lsty_arr{ i };
    scl = scl_arr( i );

    input_arr = fs.calc_forcing_term( t_arr_rhy( 1:end-1 ), W, 0, scl*eye( 2 ) );
    [ y_arr_rhy, ~, ~ ] = trans_sys.rollout( scl*q_des( :, 1 ), zeros( 2, 1), zeros( 2, 1 ), input_arr, 0, t_arr_rhy  );
    plot( a2, t_arr_rhy, y_arr_rhy( 1, : ), 'linewidth', lw, 'color', c_arr( i, : ), 'linestyle', lsty );
    plot( a4, t_arr_rhy, y_arr_rhy( 2, : ), 'linewidth', lw, 'color', c_arr( i, : ), 'linestyle', lsty );

end

ylabel( a1, 'Joint 1 (rad)', 'fontsize', 30 )
ylabel( a3, 'Joint 2 (rad)', 'fontsize', 30 )

xlabel( a3, '$t$ (s)', 'fontsize', 30 )
xlabel( a4, '$t$ (s)', 'fontsize', 30 )


set( a1, 'xlim', [ 0, 4 ], 'ylim', [-1.0, 1.0], 'fontsize', 30, 'xtick', 0:2:4 )
set( a2, 'xlim', [ 0, 4 ], 'ylim', [-1.0, 1.0], 'fontsize', 30, 'xtick', 0:2:4 )
set( a3, 'xlim', [ 0, 4 ], 'ylim', [-1.0, 1.0], 'fontsize', 30, 'xtick', 0:2:4 )
set( a4, 'xlim', [ 0, 4 ], 'ylim', [-1.0, 1.0], 'fontsize', 30, 'xtick', 0:2:4 )

fig_save( f, 'ThesisImages/images/joint_space_imit_learn_rhy_temp_scl' )

%%  -- (1C) For Both Discrete and Rhythmic Movement

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
data = load( 'learned_parameters/M2.mat' ); data = data.data;
data_ref = load( 'alphabets/M2.mat' ); 

x_real = data_ref.data.p_data( 1, : ) - data_ref.data.p_data( 1, 1 );
y_real = data_ref.data.p_data( 2, : ) - data_ref.data.p_data( 2, 1 );

N = size( data.weight, 2 );
W = data.weight;

% Append a third weight array which is zero
W( 3, 1:N ) = zeros( 1, N );
data.z0     = [ data.z0  ; 0 ];
data.goal   = [ data.goal; 0 ];

% DMP, three primitives
cs        = CanonicalSystem( 'discrete', data.tau, data.alpha_s );
trans_sys = TransformationSystem( data.alpha_z, data.beta_z, cs );
fs        = NonlinearForcingTerm( cs, N );

% Rollout with the weight array 
t0i   = 0.0;
T     = data.tau;
Nt    = 3000;
t_arr_dis = linspace( 0, T, Nt+1 );

input_arr = fs.calc_forcing_term( t_arr_dis( 1:end-1 ), W, t0i, eye( 3 ) );
[ y_arr_dis, ~, ~ ] = trans_sys.rollout( zeros( 3, 1 ), data.z0 * data.tau, data.goal, input_arr, t0i, t_arr_dis  );

f = figure( ); a = axes( 'parent', f );
hold on; 
plot( a, y_arr_dis( 1, : ), y_arr_dis( 2, : ), 'linewidth', 5, 'color', c_blue )
plot( a, x_real, y_real, 'linewidth', 8, 'linestyle', '--', 'color', 'k' )


scatter( a, y_arr_dis( 1, 1   ), y_arr_dis( 2, 1   ), 500, 'filled', 'o','markerfacecolor', [ 1,1,1], 'markeredgecolor', [0,0,0], 'linewidth', 5 )
scatter( a, y_arr_dis( 1, end ), y_arr_dis( 2, end ), 500, 'filled', 'd', 'markerfacecolor', [ 1,1,1], 'markeredgecolor', [0,0,0], 'linewidth', 5 )
axis equal
set( a, 'fontsize', 30, 'xlim', [-2,11], 'ylim', [-1,8], 'ytick', 0:2:8, 'xtick', 0:2:10 )
xlabel( a, 'X (m)', 'fontsize', 35 )
ylabel( a, 'Y (m)', 'fontsize', 35 )
title( a, 'Discrete', 'fontsize', 45 )
fig_save( f, 'ThesisImages/images/task_space_position_discrete' )


% Scale the trajectory 
scl    = [ 0.2, 0.5, 1.0, 1.5, 2.0 ];
lw_arr = [   5,   5,   8,   5,   5 ];
mk_arr = [  150, 200, 300, 350, 400 ];

f = figure(); a = axes( 'parent', f );
hold on;
for i = 1 : length( scl )
    s  = scl( i );
    lw = lw_arr( i );
    mk = mk_arr( i );

    if  abs(s-1.0) >= 1e-8
        lc = [ 0.0000, 0.0000, 0.0000 ];
    else
        lc = [ 0.0000, 0.4470, 0.7410 ];
    end

    [ y_tmp, ~, ~ ] = trans_sys.rollout( zeros( 3, 1 ), data.z0 * data.tau, s*data.goal, s*input_arr, t0i, t_arr_dis  );
    
    plot( a, y_tmp( 1, : ), y_tmp( 2, : ), 'linewidth', lw, 'color', lc )
    scatter( a, y_tmp( 1, 1   ), y_tmp( 2, 1   ), mk, 'filled', 'o','markerfacecolor', [ 1,1,1], 'markeredgecolor', [0,0,0], 'linewidth', 5 )
    scatter( a, y_tmp( 1, end ), y_tmp( 2, end ), mk, 'filled', 'd', 'markerfacecolor', [ 1,1,1], 'markeredgecolor', [0,0,0], 'linewidth', 5 )    

end
set( a, 'fontsize', 30, 'xlim', [-2,20], 'ylim', [-2, 16], 'ytick', 0:5:15, 'xtick', 0:5:20 )
set( a, 'xticklabel', {}, 'yticklabel', {}, 'zticklabel', {} )
xlabel( 'X (m)', 'fontsize', 45 )
ylabel( 'Y (m)', 'fontsize', 45 )
fig_save( f, 'ThesisImages/images/task_space_position_scale' )

% Rotate the trajectory
ang_arr = 0:60:359;
f = figure( ); a = axes( 'parent', f );
hold on
axis equal

for ang = ang_arr
    R = rotx( ang );

    if  abs( ang ) >= 1e-8
        lc = [ 0.0000, 0.0000, 0.0000 ];
        lw = 4;
    else
        lc = [ 0.0000, 0.4470, 0.7410 ];
        lw = 8;
    end

    [ y_tmp, ~, ~ ] = trans_sys.rollout( zeros( 3, 1 ), R * data.z0 * data.tau, R*data.goal, R*input_arr, t0i, t_arr_dis  );
        
    plot3( a, y_tmp( 1, : ), y_tmp( 2, : ), y_tmp( 3, : ), 'linewidth', lw, 'color', lc )
    scatter3( a, y_tmp( 1, 1   ), y_tmp( 2, 1   ), y_tmp( 3, 1 ), 500, 'filled', 'o','markerfacecolor', [ 1,1,1], 'markeredgecolor', [0,0,0], 'linewidth', 5 )
    scatter3( a, y_tmp( 1, end ), y_tmp( 2, end ), y_tmp( 3, end ), 500, 'filled', 'd', 'markerfacecolor', [ 1,1,1], 'markeredgecolor', [0,0,0], 'linewidth', 5 )    
   
end

set( a, 'view', [31.7014, 8.9169], 'fontsize',30, 'xlim', [-2, 12], 'ylim', [-8, 8], 'zlim', [-8, 8] )
set( a, 'xticklabel', {}, 'yticklabel', {}, 'zticklabel', {} )
xlabel( a, 'X (m)', 'fontsize', 45 );
ylabel( a, 'Y (m)', 'fontsize', 45 );
zlabel( a, 'Z (m)', 'fontsize', 45 );

fig_save( f, 'ThesisImages/images/task_space_position_rot' )

%%  -- (2B) For Rhythmic Movement

% Drawing a 2D Heart Symbol
% Need symbolic form
syms t

% Position
x = 16 * sin( t )^3;
y = 13*cos(t) - 5*cos(2*t) - 2*cos(3*t) - cos(4*t);

% Velocity and Acceleration
dx  = diff(  x, t ); dy  = diff(  y, t );
ddx = diff( dx, t ); ddy = diff( dy, t );

% Reformulate as functions
  x = matlabFunction(   x );   y = matlabFunction(   y );
 dx = matlabFunction(  dx );  dy = matlabFunction(  dy );
ddx = matlabFunction( ddx ); ddy = matlabFunction( ddy );

% Sample the data and use it for training
P = 100;
t_P = linspace( 0, 2*pi, P );

  x_data =   x( t_P );   y_data =   y( t_P );
 dx_data =  dx( t_P );  dy_data =  dy( t_P );
ddx_data = ddx( t_P ); ddy_data = ddy( t_P );

% Collecting the data as p_des
  p_des = [   x_data;   y_data ];
 dp_des = [  dx_data;  dy_data ];
ddp_des = [ ddx_data; ddy_data ];

% The period is set to be 2pi,hence tau = 1.0;
tau = 1;

% Parameters of DMP
% The number of Basis Function for the Nonlinear Forcing Term
% The algorithm is sensitive to the choice of N!! Care is required.
N = 25;

% Parameters of the Transformation System
alpha_z  = 50.0;
 beta_z  = 0.25 * alpha_z;

% DMP
cs        = CanonicalSystem( 'rhythmic', tau, 1.0 );
trans_sys = TransformationSystem( alpha_z, beta_z, cs );
fs        = NonlinearForcingTerm( cs, N );

% Calculating the required Nonlinear Forcing Term
% Goal location is zero 
B = trans_sys.get_desired( p_des, dp_des, ddp_des, mean( p_des,2 ) );

% The A matrix 
A = zeros( N, P );

% Interating along the sample points
for i = 1:P 
    t = t_P( i );
    A( :, i ) = fs.calc_multiple_ith( t, 1:N )/ fs.calc_whole_at_t( t );
end

% Linear Least-square
W = B * A' * ( A * A' )^-1;

% Rollout with the weight array 
t0i   = 0.0;
T     = tau*9;
Nt    = 3000;
t_arr_dis = linspace( 0, T, Nt+1 );

scl = 1.0;
input_arr = fs.calc_forcing_term( t_arr_dis( 1:end-1 ), scl*W, t0i, eye( 2 ) );
[ y_arr_dis, ~, ~ ] = trans_sys.rollout( scl*p_des( :, 1 ), dp_des( :, 1 )*tau, scl*mean( p_des,2 ), input_arr, t0i, t_arr_dis  );

f = figure( ); a = axes( 'parent', f );
hold on; 
plot( a, y_arr_dis( 1, : ), y_arr_dis( 2, : ), 'linewidth', 5, 'color', c_blue )
plot( a, x_data, y_data, 'linewidth', 8, 'color', 'k', 'linestyle', '--' );
scatter( a, y_arr_dis( 1, 1   ), y_arr_dis( 2, 1   ), 500, 'filled', 'o','markerfacecolor', [ 1,1,1], 'markeredgecolor', [0,0,0], 'linewidth', 5 )

axis equal
set( a, 'fontsize', 30, 'xlim', [-20,20], 'ylim', [-20,15] )
xlabel( a, 'X (m)', 'fontsize', 35 )
ylabel( a, 'Y (m)', 'fontsize', 35 )
title( a, 'Rhythmic', 'fontsize', 45 )

fig_save( f, 'ThesisImages/images/task_space_position_rhythmic' )

%% =======================================================
%% (2-) Imitation Learning for Task-space, Orientation, Section 4.4.3
%%  -- (2A) For Discrete Movement

% Generate an example orientation trajectory 
Ri = eye( 3 );
Rg = rotx( 30 ) * roty( 20 ) * rotz( 50 );

% Get the axis (or diff between two SO3 matrices)
w = LogSO3( Ri' * Rg );
theta = norm( w );
w_hat = w/theta;

% We use the symbolic form to get the angular velocity 
syms t positive

% Get the symbolic form for the trajectory
D = 3.0;
traj = 10*(t/D)^3 - 15*(t/D)^4 + 6*(t/D)^5;
R_sym = Ri * ( eye( 3 ) + sin( theta*traj ) * w_hat + (1 - cos( theta*traj ) ) * w_hat^2 ); 

% Also get the angular velocity in symbolic form
ang_vel_sym = diff( R_sym, t ) * R_sym.';

% Change these to MATLAB_functions
R_func = matlabFunction( R_sym );
ang_vel_func = matlabFunction( ang_vel_sym );

% Regenerate the trajectories 
P = 500;
t_des = linspace( 0, D, P );

R_des   = zeros( 3, 3, P );
ang_des = zeros( 3, 3, P );

for i = 1 : P
   t = t_des( i ); 
   R_des(   :, :, i ) = R_func( t );
   ang_des( :, :, i ) = ang_vel_func( t );
end

% Calculate the error vector and its derivative
  err = zeros( 3, P );
 derr = zeros( 3, P );
dderr = zeros( 3, P );

for i = 1 : P
    R   = R_des( :, :, i );
    ang = ang_des( :, :, i );
    R0g = R' * Rg;
   err( :, i ) = so3_to_R3( LogSO3( R0g ) );
   
   dR = -R' * ang * Rg;
    
   % Get theta
   th = norm( LogSO3( R0g ) );
   
   if th >= 1e-7
        tmp1 = ( ( th * cos( th ) - sin( th ) ) / ( 4 * sin( th )^3 ) ) * trace( dR ) * ( R0g - R0g' );
        tmp2 = -th/( 2 * sin( th ) ) * ( dR - dR' );
   else
      tmp1 = zeros( 3, 3 );
      tmp2 = zeros( 3, 3 );
   end
   
   derr( :, i ) = so3_to_R3( tmp1 + tmp2 );
   
end

% For dderr, do numerical differentiation
for i = 1 : 3
    dderr( i, : ) =  data_diff( derr( i, : ) );
    dderr( i, : ) = smoothdata( dderr( i, : ), "gaussian", 50 );
end


% Parameters of the Transformation System
tau       = D;
alpha_s   = 1.0;
alpha_z   = 500.0;
beta_z    = 0.25 * alpha_z;
N         = 30;

% DMP, three primitives
cs        = CanonicalSystem( 'discrete', tau, alpha_s );
trans_sys = TransformationSystem( alpha_z, beta_z, cs );
fs        = NonlinearForcingTerm( cs, N );

% Scaling vector 
Se = diag( so3_to_R3( LogSO3( Ri' * Rg ) ) );

% Calculating the required Nonlinear Forcing Term
B = inv( Se ) * trans_sys.get_desired( err, derr, dderr, zeros( 3, 1 ) );


% The A matrix 
A = zeros( N, P );

% Interating along the sample points
for i = 1:P 
    t = t_des( i );
    A( :, i ) = fs.calc_multiple_ith( t, 1:N )/ fs.calc_whole_at_t( t ) * cs.calc( t );
end

% Linear Least-square
W = B * A' * ( A * A' )^-1;

% Rollout with the weight array 
t0i   = 0.0;
T     = D*2;
Nt    = 3000;
t_arr_dis = linspace( 0, T, Nt+1 );


% Scaling vector 
Se2 = diag( so3_to_R3( LogSO3( Ri' * Rg ) ) );

input_arr = fs.calc_forcing_term( t_arr_dis( 1:end-1 ), W, t0i, Se2 );
[ err_roll, ~, ~ ] = trans_sys.rollout( err( :, 1 ), derr( :, 1 )*tau, zeros( 3, 1 ), input_arr, t0i, t_arr_dis  );

f= figure(); a = axes( 'parent', f);
hold on
plot( a, t_arr_dis, err_roll,'linewidth', 5, 'color', c_blue )
plot( a, t_des, err,'linewidth', 8, 'color', 'k', 'linestyle', '--' )
