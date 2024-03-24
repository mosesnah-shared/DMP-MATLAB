% [Title]     Plotting the Nonlinear Forcing Term
% [Author]    Moses Chong-ook Nah
% [Email]     mosesnah@mit.edu
% [Update]    At 2023.11.12

%% (--) Initialization
clear; close all; clc;
fig_config( 'fontSize', 20, 'markerSize', 10 )

%% =======================================================
%% (1-) Nonlinear Forcing Term
%%  -- (1A) Nonlinear Forcing Term for Discrete Movement

dt      = 1e-3;     % Time-step
T       =  3.0;     % Total Time
tau     =    T;     % Time Constant
alpha_s =  1.0;     % The Gain of the Canonical System
N       =   30;     % Number of Basis Functions
t_arr   = 0:dt:T;   % Time array

% Canonical System and Nonlinear Forcing Term
type = 'discrete'; 
cs = CanonicalSystem( type, tau, alpha_s ); 
fs = NonlinearForcingTerm( cs, N );

a1 = subplot( 2, 1, 1 );
a2 = subplot( 2, 1, 2 );
hold( a1, 'on' );
hold( a2, 'on' );

plot( a1,          t_arr  , fs.calc_multiple_ith( t_arr, 1:N ), 'linewidth', 3 );
plot( a2, cs.calc( t_arr ), fs.calc_multiple_ith( t_arr, 1:N ), 'linewidth', 3 );    

% Overlap of the whole activation
plot( a1, t_arr, fs.calc_whole_at_t( t_arr ), 'color', 'k' )

set( a1, 'xticklabel', { }, 'yticklabel', { }, 'xlim', [ 0, T ] )
set( a2, 'xticklabel', { }, 'yticklabel', { }, 'xlim', [ 0, 1 ] )


%%  -- (1B) Nonlinear Forcing Term for Rhythmic Movement

% Canonical System and Nonlinear Forcing Term
type = 'rhythmic'; 


t_arr   = 0:1e-3:8.0;
Tp      =   4.0;         % Period of the Rhythmic movement
tau     =   Tp/(2*pi);   % Time Constant, which is Tp/2pi
alpha_s = 1.0;
N = 5;

cs = CanonicalSystem( type, tau, alpha_s ); 
fs = NonlinearForcingTerm( cs, N );

a1 = subplot( 2, 1, 1 );
a2 = subplot( 2, 1, 2 );
hold( a1, 'on' );
hold( a2, 'on' );

plot( a1,          t_arr  , fs.calc_multiple_ith( t_arr, 1:N ), 'linewidth', 3 );
plot( a2, cs.calc( t_arr ), fs.calc_multiple_ith( t_arr, 1:N ), 'linewidth', 3 );    

% Overlap of the whole activation
% plot( a1, t_arr, fs.calc_whole_at_t( t_arr ), 'color', 'k' )
plot( a2, cs.calc( t_arr ), fs.calc_whole_at_t( t_arr ), 'color', 'k', 'linewidth', 5 )
