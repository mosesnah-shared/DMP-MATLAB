%% [Example Script] Plotting the Nonlinear Forcing Terms
% [Author] Moses Chong-ook Nah
% [Date]   2023.08.15

%% [Initialization] 
clear; close all; clc;
fig_config( 'fontSize', 20, 'markerSize', 10 )

%% [Nonlinear Forcing Term] Discrete
tau     = 1.0;
alpha_s = 1.0;
cs_discrete = CanonicalSystem( 'discrete', tau, alpha_s ); 

N = 10;
g = 1.0; y0 = 0.0; r = 1.0; % 3rd argument is ignored for discrete movement

fs = NonlinearForcingTerm( cs_discrete, N );
t = 0:0.001:2.0;

f = figure( ); a = axes( 'parent', f );
hold on
for i = 1:N
    plot( t, fs.calc_ith( t, i ) );
%     plot( cs_discrete.calc( t ), fs.calc_ith( t, i ) );
end

%% [Nonlinear Forcing Term] Rhythmic 
% Will develop once necessary
