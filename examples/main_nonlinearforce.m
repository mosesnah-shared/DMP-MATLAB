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

w_arr = 2 *rand( 1, N );

f = figure( ); a = axes( 'parent', f );
hold on
for i = 1:N
    plot( t, fs.calc_ith( t, i ), 'color', [0 0.4470 0.7410] );
%     plot( cs_discrete.calc( t ), fs.calc_ith( t, i ), 'color', [0 0.4470 0.7410] );
end

set(gca, 'xlim', [0,1.2])
set( gca, 'xticklabel', {}, 'yticklabel', {} )


%% 
% Also plot the weighted version
tmp = zeros( 1, length( t ) );
for i = 1 : N
    tmp = tmp + fs.calc_ith( t, i ) .* w_arr( i ) .* cs_discrete.calc( t ); 
end

plot( t, tmp, 'color', 'k' )

set(gca, 'xlim', [0,1.3])

%% [Nonlinear Forcing Term] Rhythmic 
% Will develop once necessary


%% [Nonlinear Forcing Term] Discrete
tau     = 1.0;
alpha_s = 1.0;
cs_rhythmic = CanonicalSystem( 'rhythmic', tau, alpha_s ); 

N = 10;
g = 1.0; y0 = 0.0; r = 1.0; % 3rd argument is ignored for discrete movement

fs = NonlinearForcingTerm( cs_rhythmic, N );
t = 0:0.001:2*pi;

w_arr = 2 *rand( 1, N );

f = figure( ); a = axes( 'parent', f );
hold on
for i = 1:N
    plot( t, fs.calc_ith( t, i ), 'color', [0 0.4470 0.7410] );
%     plot( cs_discrete.calc( t ), fs.calc_ith( t, i ), 'color', [0 0.4470 0.7410] );
end

set(gca, 'xlim', [0,2*pi])
set( gca, 'xticklabel', {}, 'yticklabel', {} )


