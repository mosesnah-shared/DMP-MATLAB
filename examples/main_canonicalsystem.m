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
for i = 1:4
    cs_discrete.alpha_s = i;
    plot( t, cs_discrete.calc( t ) );
end

%% [Canonical System] Rhytmic
t    = 0:0.01:30.0;
tau  = 1.0;

cs_rhythmic = CanonicalSystem( 'rhythmic', 1.0, 1.0 ); % 3rd argument ignored

f = figure( ); a = axes( 'parent', f );
hold on
cs_rhythmic.tau = 0.5;
plot( t, cs_rhythmic.calc( t ) );

