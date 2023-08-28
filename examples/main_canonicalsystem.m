%% ==================================================================
%% [Title] Example Script - Plotting the Canonical System
% Author: Moses Chong-ook Nah
%  Email: mosesnah@mit.edu
%   Date: 2023.08.15
%% ==================================================================

%% [0A] Initialization
clear; close all; clc;
fig_config( 'fontSize', 20, 'markerSize', 10 )


%% [1-] Canonical System 
%% ==== [1A] Discrete Movement

dt     = 1e-3;             % Time-step
T      =  1.0;             % Total Time
tau    =    T;             % Time Constant
alphas =  1.0;             % The gain for the canonical system

% The Time array
t_arr  = 0:dt:T;

% Construct the discrete Canonical System
cs_discrete = CanonicalSystem( 'discrete', tau, alphas );

% Plotting the time vs. canonical system with respect to alpha_s
f = figure( ); a = axes( 'parent', f );
hold( a, 'on' )

% Iterating over the alphas values
for as = [ 0.1, 0.2, 0.5, 1.0, 2.0 ]
    cs_discrete.alpha_s = as;
    plot( a, t_arr, cs_discrete.calc( t_arr ), 'linewidth', 5 );
    set( a, 'xticklabel', {}, 'yticklabel', {}, 'xlim', [0, T] )
end

title( a, 'Discrete Canonical System', 'fontsize', 30)

%% ==== [1B] Rhythmic Movement

dt     =  1e-3;             % Time-step
T      =  20.0;             % Total Time
tau    =  4.0/(2*pi);       % Time Constant, which is the period/2pi

% The Time array
t_arr  = 0:dt:T;

% Third argument is ignored since 'rhythmic' movement 
% does not need that parameter
cs_rhythmic = CanonicalSystem( 'rhythmic', tau, alphas ); % 

f = figure( ); a = axes( 'parent', f );
hold on
plot( a, t_arr, cs_rhythmic.calc( t_arr ), 'linewidth', 5 );
set( a, 'xticklabel', {}, 'yticklabel', {}, 'xlim', [0,T] )
title( a, 'Rhythmic Canonical System', 'fontsize', 30)

