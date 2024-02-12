% [Title]     Generating Images for Contraction Theory Paper
% [Author]    Moses Chong-ook Nah
% [Email]     mosesnah@mit.edu
% [Update]    At 2024.02.12

%% (--) Initialization
clear; close all; clc;
fig_config( 'fontSize', 20, 'markerSize', 10 )

%% (1A) Dynamic Movement Primitives for multiple movements

% For discrete movement, load the data
tmp  = load( './learned_parameters/discrete/A.mat' );
data = tmp.data;

% Getting the number of basis functions from the nonlinear forcing term
[ ~, N ] = size( data.weight );

cs1        = CanonicalSystem( 'discrete', data.tau, data.alpha_s );
fs1        = NonlinearForcingTerm( cs1, N );
trans_sys1 = TransformationSystem( data.alpha_z, data.beta_z, cs1 );

% The parameters for forward simulation
t0i   = 0.0;
T     = 16;
dt    = 1e-2;
t_arr = 0:dt:T;
Nt    = length( t_arr );

% Calculate the nonlinear forcing term for discrete movement
input_arr = fs1.calc_forcing_term( t_arr( 1:end-1 ), data.weight, t0i, eye( 2 ), 'trimmed' );

% Rollout, note that initial condition is set to be zeros. 
[ y_arr, z_arr, dy_arr ] = trans_sys1.rollout( zeros( 2, 1 ), zeros( 2, 1 ), data.goal, input_arr, t0i, t_arr  );  

f = figure( ); a = axes( 'parent', f );
a1 = subplot( 3, 2, 1 );
plot( a1, t_arr, cs1.calc( t_arr ), 'linewidth', 4, 'color',  [0 0.4470 0.7410] );
set( a1, 'xticklabel', {}, 'yticklabel', {}, 'xlim', [0,16] )
xlabel( '$t$ (s)', 'fontsize', 40 )
ylabel( '$s_d$', 'fontsize', 40 )

a3 = subplot( 3, 2, 3 );
hold on
plot( a3, t_arr( 1:end-1 ), input_arr( 1, :), 'linewidth', 4, 'linestyle', '-' , 'color', [0 0.4470 0.7410]	);
plot( a3, t_arr( 1:end-1 ), input_arr( 2, :), 'linewidth', 4, 'linestyle', '-' , 'color', [0 0.4470 0.7410]	);
set( a3, 'xticklabel', {}, 'yticklabel', {}, 'xlim', [0,6.0], 'ylim', 0.3*[-1e+6,1e+6] )
xlabel( '$t$ (s)', 'fontsize', 40 )
ylabel( '$\mathbf{F}_d$', 'fontsize', 40 )

a5 = subplot( 3, 2, 5 );
hold on
plot( a5, y_arr( 1, : ), y_arr( 2, : ), 'linewidth', 4, 'color', [0 0.4470 0.7410] );
scatter( a5, y_arr( 1,   1 ), y_arr( 2,   1 ), 500, 'filled',  'o', 'markerfacecolor', 'w', 'markeredgecolor', [0 0.4470 0.7410], 'linewidth', 4 )
scatter( a5, y_arr( 1, end ), y_arr( 2, end ), 500, 'filled',  'd', 'markerfacecolor', 'w', 'markeredgecolor', [0 0.4470 0.7410], 'linewidth', 4 )
set( a5, 'xticklabel', {}, 'yticklabel', {}, 'xlim', [-3, 9.5], 'ylim', [-2, 10.5] )
axis equal
xlabel( '$X$ (m)', 'fontsize', 40 )
ylabel( '$Y$ (m)', 'fontsize', 40 )

% For rhythmic movement, load the data
tmp  = load( './learned_parameters/rhythmic/heart.mat' );
data = tmp.data;

% Getting the number of basis functions from the nonlinear forcing term
[ ~, N ] = size( data.weight );

cs2        = CanonicalSystem( 'rhythmic', data.tau, 1.0 );
fs2        = NonlinearForcingTerm( cs2, N );
trans_sys2 = TransformationSystem( data.alpha_z, data.beta_z, cs2 );

% The parameters for forward simulation
t0i   = 0.0;
T     = 3;
dt    = 1e-4;
t_arr = 0:dt:T;
Nt    = length( t_arr );

input_arr2 = fs2.calc_forcing_term( t_arr( 1:end-1 ), data.weight , t0i, eye( 2 ) );
[ y_arr2, ~, ~ ] = trans_sys2.rollout( data.p_init, zeros( 2, 1 ), [ 0; 0.8 ], input_arr2, t0i, t_arr  );

a2 = subplot( 3, 2, 2 );
plot( t_arr, cs2.calc( t_arr ), 'linewidth', 4, 'color', [0.8500 0.3250 0.0980] );
set( a2, 'xticklabel', {}, 'yticklabel', {}, 'xlim', [0,T] )
xlabel( '$t$ (s)', 'fontsize', 40 )
ylabel( '$s_r$', 'fontsize', 40 )

a4 = subplot( 3, 2, 4 );
hold on
plot( a4, t_arr( 1:end-1 ), input_arr2( 1, :), 'linewidth', 4, 'linestyle', '-' , 'color', [0.8500 0.3250 0.0980] );
plot( a4, t_arr( 1:end-1 ), input_arr2( 2, :), 'linewidth', 4, 'linestyle', '-' , 'color', [0.8500 0.3250 0.0980] );
set( a4, 'xticklabel', {}, 'yticklabel', {}, 'xlim', [0,T], 'ylim', [-6*10^5,6*10^5] )
xlabel( '$t$ (s)', 'fontsize', 40 )
ylabel( '$\mathbf{F}_r$', 'fontsize', 40 )

a6 = subplot( 3, 2, 6 );

plot( a6, y_arr2( 1, 1:end-1 ), y_arr2( 2, 1:end-1 ), 'linewidth', 4, 'color', [0.8500 0.3250 0.0980] );
set( a6, 'xticklabel', {}, 'yticklabel', {}, 'xlim', [-20,20], 'ylim',[-20,20])
axis equal
xlabel( '$X$ (m)', 'fontsize', 40 )
ylabel( '$Y$ (m)', 'fontsize', 40 )

export_fig ./Images/DMP_images/dis_and_rhy_DMP.pdf -transparent
