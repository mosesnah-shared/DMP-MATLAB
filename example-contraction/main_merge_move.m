% [Title]     Contraction Theory, Merging Movements
% [Author]    Moses Chong-ook Nah
% [Email]     mosesnah@mit.edu
% [Update]    At 2024.01.31

%% (--) Initialization
clear; close all; clc;
fig_config( 'fontSize', 20, 'markerSize', 10 )

%% =======================================================
%% (1-) Merging Two Discrete Movements
%%  -- (1A) Calling the Weight Values and check the output
close all;

% Load two mat functions, discrete movements
tmp1 = load( './learned_parameters/discrete/M.mat' ); dis1 = tmp1.data;
tmp2 = load( './learned_parameters/discrete/A.mat' ); dis2 = tmp2.data;

% Get the degrees of freedom and the number of weights 
n = size( dis1.weight, 1 );
N = size( dis1.weight, 2 );

% The alpha, beta values must match
% The Three elements of DMP
cs        = CanonicalSystem( 'discrete', dis1.tau, dis1.alpha_s );
fs        = NonlinearForcingTerm( cs, N );
trans_sys = TransformationSystem( dis1.alpha_z, dis1.beta_z, cs );

% Rollout for each discrete movement
t0i   =    1.0;           % The initial time of the movement rollout
T     =    5.0;           % The   whole time of the simulation 
dt    =   1e-3;           % Time-step  for the simulation
t_arr = 0:dt:T;           % Time array for the simulation

F_arr1 = fs.calc_forcing_term( t_arr( 1:end-1 ), dis1.weight, t0i, eye( 2 ) );
F_arr2 = fs.calc_forcing_term( t_arr( 1:end-1 ), dis2.weight, t0i, eye( 2 ) );

[ y_arr1, ~, ~] = trans_sys.rollout( zeros( n, 1 ), zeros( n, 1 ), dis1.goal, F_arr1, t0i, t_arr  );    
[ y_arr2, ~, ~] = trans_sys.rollout( zeros( n, 1 ), zeros( n, 1 ), dis2.goal, F_arr2, t0i, t_arr  );    
 
f = figure( ); a = axes( 'parent', f );
hold on; axis equal
plot( a, y_arr1( 1, : ), y_arr1( 2, : ), 'linewidth', 6 )
plot( a, y_arr2( 1, : ), y_arr2( 2, : ), 'linewidth', 6 )

plot( a, 0.5*y_arr1( 1, : )+0.5*y_arr2( 1, : ), 0.5*y_arr1( 2, : )+0.5*y_arr2( 2, : ), 'linewidth', 6 )
plot( a, y_arr2( 1, : ), y_arr2( 2, : ), 'linewidth', 6 )

%%  -- (1B) Draw the Results
close all;
f = figure( ); a = axes( 'parent', f );
hold on; axis equal

gain_arr = 0:0.1:1;
Ntmp = length( gain_arr );
off = 10.0;
for i = 1: Ntmp
    gain = gain_arr( i );

    F_arr =    F_arr1 * gain +    F_arr2 * ( 1 - gain );
    goal  = dis1.goal * gain + dis2.goal * ( 1 - gain );

    [ y_arr_comb, ~, ~] = trans_sys.rollout( zeros( n, 1 ), zeros( n, 1 ), goal, F_arr, t0i, t_arr  );    
    
    plot( a, (i-1)*off+y_arr_comb( 1, : ), y_arr_comb( 2, : ), 'linewidth', 6 )
end

%% =======================================================
%% (2-) Merging Discrete and Rhythmic Movements
%%  -- (2A) Calling the Weight Values and check the output
close all;

% Load two mat functions, discrete movements
tmp1 = load( './learned_parameters/discrete/M.mat' );      dis = tmp1.data;
tmp2 = load( './learned_parameters/rhythmic/circle.mat' ); rhy = tmp2.data;

% Get the degrees of freedom and the number of weights 
n = size( dis.weight, 1 );
N = size( dis.weight, 2 );

% Create DMP for discrete
cs_d        = CanonicalSystem( 'discrete', dis.tau*2, dis.alpha_s );
fs_d        = NonlinearForcingTerm( cs_d, N );
trans_sys_d = TransformationSystem( dis.alpha_z, dis.beta_z, cs_d );

% Create DMP for rhythmic
cs_r        = CanonicalSystem( 'rhythmic', rhy.tau, 1.0 );
fs_r        = NonlinearForcingTerm( cs_r, N );
trans_sys_r = TransformationSystem( rhy.alpha_z, rhy.beta_z, cs_r );

% Rollout for each discrete movement
t0i   =    0.0;           % The initial time of the movement rollout
T     =   10.0;           % The   whole time of the simulation 
dt    =   1e-3;           % Time-step  for the simulation
t_arr = 0:dt:T;           % Time array for the simulation

scl_r = 8.0;

F_arr_d = fs_d.calc_forcing_term( t_arr( 1:end-1 ), dis.weight, t0i, eye( 2 ) );
F_arr_r = scl_r*fs_r.calc_forcing_term( t_arr( 1:end-1 ), rhy.weight, t0i, eye( 2 ) );

[ y_arr_d, ~, ~] = trans_sys_d.rollout( zeros( n, 1 ), zeros( n, 1 ), dis.goal, F_arr_d, t0i, t_arr  );    
[ y_arr_r, ~, ~] = trans_sys_r.rollout( scl_r*rhy.p_init, scl_r*rhy.tau*rhy.dp_init, zeros( 2, 1), F_arr_r, t0i, t_arr  );    

f = figure( ); a = axes( 'parent', f );
hold on; axis equal
plot( a, y_arr_d( 1, : ), y_arr_d( 2, : ), 'linewidth', 6 )
plot( a, y_arr_r( 1, : ), y_arr_r( 2, : ), 'linewidth', 6 )

%%  -- (2B) Draw the Results
close all;
f = figure( ); a = axes( 'parent', f );
hold on; axis equal

gain_arr = 0:0.1:1;
Ntmp = length( gain_arr );
off = 10.0;
for i = 1: Ntmp
    gain = gain_arr( i );
    F_arr =  F_arr_d * gain + F_arr_r * ( 1 - gain );
    goal  =  dis.goal * gain;

    [ y_arr_comb, ~, ~] = trans_sys_d.rollout( zeros( n, 1 ), zeros( n, 1 ), goal, F_arr, t0i, t_arr  );    
    
    plot( a, (i-1)*off+y_arr_comb( 1, : ), y_arr_comb( 2, : ), 'linewidth', 6 )
end

%% =======================================================
%% (3-) Merging Two Rhythmic Movements
%%  -- (3A) Calling the Weight Values and check the output
close all;

% Load two mat functions, discrete movements
tmp1 = load( './learned_parameters/rhythmic/heart.mat' );  rhy1 = tmp1.data;
tmp2 = load( './learned_parameters/rhythmic/circle.mat' ); rhy2 = tmp2.data;

% Get the degrees of freedom and the number of weights 
n = size( rhy1.weight, 1 );
N = size( rhy1.weight, 2 );

% Create DMP for rhythmic
cs_r        = CanonicalSystem( 'rhythmic', rhy1.tau, 1.0 );
fs_r        = NonlinearForcingTerm( cs_r, N );
trans_sys_r = TransformationSystem( rhy1.alpha_z, rhy1.beta_z, cs_r );

% Rollout for each discrete movement
T     =    5.0;           % The   whole time of the simulation 
dt    =   1e-3;           % Time-step  for the simulation
t_arr = 0:dt:T;           % Time array for the simulation

scl1 = 1;
scl2 = 20;

F_arr1 = fs_r.calc_forcing_term( t_arr( 1:end-1 ), rhy1.weight, 0, scl1 * eye( 2 ) );
F_arr2 = fs_r.calc_forcing_term( t_arr( 1:end-1 ), rhy2.weight, 0, scl2 * eye( 2 ) );

[ y_arr1, ~, ~] = trans_sys.rollout( scl1 * rhy1.p_init, scl1 * rhy1.dp_init*rhy1.tau, zeros( 2, 1 ), F_arr1, 0, t_arr  );    
[ y_arr2, ~, ~] = trans_sys.rollout( scl2 * rhy2.p_init, scl2 * rhy2.dp_init*rhy2.tau, zeros( 2, 1 ), F_arr2, 0, t_arr  );    
 
f = figure( ); a = axes( 'parent', f );
hold on; axis equal
plot( a, y_arr1( 1, : ), y_arr1( 2, : ), 'linewidth', 6 )
plot( a, y_arr2( 1, : ), y_arr2( 2, : ), 'linewidth', 6 )

%%  -- (2B) Draw the Results

close all;
f = figure( ); a = axes( 'parent', f );
hold on; axis equal

gain_arr = 0:0.1:1;
Ntmp = length( gain_arr );
off = 10.0;

for i = 1: Ntmp
    gain = gain_arr( i );

    F_arr =    F_arr1 * gain +    F_arr2 * ( 1 - gain );
     init = scl1 *  rhy1.p_init * gain + scl2 * rhy2.p_init * (1 - gain )
    dinit = scl1 * rhy1.dp_init * gain + scl2 * rhy2.dp_init * (1 - gain )

    [ y_arr_comb, ~, ~] = trans_sys.rollout( init, dinit, zeros( 2, 1), F_arr, 0, t_arr  );    
    
    plot( a, (i-1)*off+y_arr_comb( 1, : ), y_arr_comb( 2, : ), 'linewidth', 6 )
end
