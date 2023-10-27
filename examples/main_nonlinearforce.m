%% ==================================================================
%% [Title] Example Script - Plotting the Nonlinear Forcing Terms
% Author: Moses Chong-ook Nah
%  Email: mosesnah@mit.edu
%   Date: 2023.08.15
%% ==================================================================

%% [0A] Initialization
clear; close all; clc;
fig_config( 'fontSize', 20, 'markerSize', 10 )

%% [1-] Nonlinear Forcing Term
%% ==== [1A] Discrete/Rhythmic Movement, Individual and Whole Activation Functions

dt      = 1e-3;             % Time-step
T       =  3.0;             % Total Time
tau     =    T;             % Time Constant
alpha_s =  1.0;             % The gain for the canonical system
N       =   10;             % Number of Basis Functions

% Time Array
t_arr = 0:dt:T;

% Canonical System and Nonlinear Forcing Term
type = 'discrete'; %'rhythmic', If rhythmic, must adapt variable tau 
cs = CanonicalSystem( type, tau, alpha_s ); 
fs = NonlinearForcingTerm( cs, N );

a1 = subplot( 2, 1, 1 );
a2 = subplot( 2, 1, 2 );
hold( a1, 'on' );
hold( a2, 'on' );

plot( a1,          t_arr  , fs.calc_ith_arr( t_arr, 1:N ), 'color', 'r' );
plot( a2, cs.calc( t_arr ), fs.calc_ith_arr( t_arr, 1:N ), 'color', 'b' );    

% Overlap the whole activation
plot( a1, t_arr, fs.calc_whole_at_t( t_arr ), 'color', 'k' )

set( a1, 'xticklabel', {}, 'yticklabel', {}, 'xlim', [ 0, T ] )
set( a2, 'xticklabel', {}, 'yticklabel', {}, 'xlim', [ 0, 1 ] )

