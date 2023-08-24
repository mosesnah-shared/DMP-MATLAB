%% [Example Script] Plotting the Canonical System
% [Author] Moses Chong-ook Nah
% [Date]   2023.08.15

%% [Initialization] 
clear; close all; clc;
fig_config( 'fontSize', 20, 'markerSize', 10 )

%% [Canonical System] Discrete
t   = 0:0.001:1.0;
tau = 1.0;
cs_discrete = CanonicalSystem( 'discrete', tau, 1.0 );

f = figure( ); a = axes( 'parent', f );
hold on
cs_discrete.alpha_s = i;
plot( t, cs_discrete.calc( t ), 'color', [0.8500 0.3250 0.0980] );
set( gca, 'xticklabel', {}, 'yticklabel', {} )

%% [Canonical System] Rhytmic
t    = 0:0.01:44.0;
tau  = 1.0;

cs_rhythmic = CanonicalSystem( 'rhythmic', 1.0, 1.0 ); % 3rd argument ignored

f = figure( ); a = axes( 'parent', f );
hold on
cs_rhythmic.tau = 0.5*pi;
plot( t, cs_rhythmic.calc( t ),  'color', [0.8500 0.3250 0.0980] );
set( gca, 'xticklabel', {}, 'yticklabel', {}, 'xlim', [0,44] )

