% [Title]     Plotting the Canonical System
% [Author]    Moses Chong-ook Nah
% [Email]     mosesnah@mit.edu
% [Update]    At 2023.11.12

%% (--) Initialization
clear; close all; clc;
fig_config( 'fontSize', 20, 'markerSize', 10 )

%% =======================================================
%% (1-) Plotting the response of the Canonical System
%%  -- (1A) Canonical System for Discrete Movement

dt     = 1e-3;     % Time-step
T      =  1.0;     % Total Time
tau    =    T;     % Time Constant of the Canonical System
alphas =  1.0;     % The gain for the Canonical System
t_arr  = 0:dt:T;   % Time array

% Construct the discrete Canonical System
cs_discrete = CanonicalSystem( 'discrete', tau, alphas );

% Plotting the time vs. canonical system with respect to alpha_s
f = figure( ); a = axes( 'parent', f );
hold( a, 'on' )

% Iterating over the alphas values
for as = 1:2:20
    cs_discrete.alpha_s = as;
    plot( a, t_arr, cs_discrete.calc( t_arr ), 'linewidth', 5 );
    set( a, 'xticklabel', { }, 'yticklabel', { }, 'xlim', [0, T] )
end

title( a, 'Discrete Canonical System', 'fontsize', 30 )

%%  -- (1B) Canonical System for Rhythmic Movement

dt     =  1e-3;         % Time-step
T      =  20.0;         % Total Time
Tp     =   4.0;         % Period of the Rhythmic movement
tau    =   Tp/(2*pi);   % Time Constant, which is the Tp/2pi
t_arr  = 0:dt:T;        % The Time array

% Third argument is ignored since 'rhythmic' movement 
% does not need that parameter
cs_rhythmic = CanonicalSystem( 'rhythmic', tau, 1.0 ); 

f = figure( ); a = axes( 'parent', f );
hold on
plot(  a, t_arr, cs_rhythmic.calc( t_arr ), 'linewidth', 5 );
set(   a, 'xticklabel', { }, 'yticklabel', { }, 'xlim', [0,T] )
title( a, 'Rhythmic Canonical System', 'fontsize', 30)

