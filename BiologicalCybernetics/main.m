% [Title]     Generating Images for Contraction Theory Paper
% [Author]    Moses Chong-ook Nah
% [Email]     mosesnah@mit.edu
% [Update]    At 2024.02.12

%% (--) Initialization
clear; close all; clc;
fig_config( 'fontSize', 20, 'markerSize', 10 )

c_blue   = [0.0000 0.4470 0.7410];
c_orange = [0.8500 0.3250 0.0980];

%% (1A) [Figure 1] Discrete and Rhythmic DMPs

close all; clc; 

% ==================================== %
% Discrete DMP
% ==================================== %
tmp  = load( '../learned_parameters/discrete/A_loose.mat' ); 
data_discrete = tmp.data;

% The number of basis functions N
W_dis = data_discrete.weight; 
g_dd  = data_discrete.goal;
[ ~, N ] = size( W_dis );

% The gains for the transformation system
as = data_discrete.alpha_s;
az = data_discrete.alpha_z; 
bz = data_discrete.beta_z;

% Duration of the discrete movement, demonstrated trajectory
tau_dd = data_discrete.tau;  

% The three elements of discrete DMP, subscript d is used
sd   = CanonicalSystem( 'discrete', tau_dd, as );
fd   = NonlinearForcingTerm( sd, N );
ts_d = TransformationSystem( az, bz, sd );

% The Parameters for Forward Simulation
t0    = 0.0;
T     = 16;
dt    = 1e-2;
t_arr = 0:dt:T;

% Calculate the nonlinear forcing term for discrete movement and rollout
input_arr = fd.calc_forcing_term( t_arr( 1:end-1 ), W_dis, t0, eye( 2 ) );

% The primitive, including the goal location
prim_dis  = input_arr + az * bz * g_dd;
[ y_arr, ~, ~ ] = ts_d.rollout( zeros( 2, 1 ), zeros( 2, 1 ), zeros( 2, 1 ), prim_dis, t0, t_arr  );  

% Drawing the plots
f = figure( ); 
a1 = subplot( 3, 2, 1 );
plot( a1, t_arr, sd.calc( t_arr ), 'linewidth', 4, 'color',  c_blue );
set(  a1,'xlim', [0,6], 'fontsize', 30, 'xtick', [0:2:6], 'ytick', [0,1],'ylim', [0, 1] )
ylabel( '$s_d(t)$', 'fontsize', 40 )

a3 = subplot( 3, 2, 3 );
hold on
plot( a3, t_arr( 1:end-1 ), input_arr( 1, : ), 'linewidth', 4, 'linestyle', '-' , 'color', c_blue );
plot( a3, t_arr( 1:end-1 ), input_arr( 2, : ), 'linewidth', 4, 'linestyle', '-' , 'color', c_blue );
set(  a3, 'yticklabel', {}, 'xlim', [0,6.0], 'fontsize', 30, 'xtick', [0:2:6] )
xlabel( '$t$', 'fontsize', 40 )
ylabel( '$\mathbf{f}_d(t)$', 'fontsize', 40 )

a5 = subplot( 3, 2, 5 );
hold on
plot( a5, y_arr( 1, : ), y_arr( 2, : ), 'linewidth', 4, 'color', c_blue );
scatter( a5, y_arr( 1,   1 ), y_arr( 2,   1 ), 500, 'filled',  'o', 'markerfacecolor', 'w', 'markeredgecolor', c_blue, 'linewidth', 4 )
scatter( a5, y_arr( 1, end ), y_arr( 2, end ), 500, 'filled',  'd', 'markerfacecolor', 'w', 'markeredgecolor', c_blue, 'linewidth', 4 )
set( a5, 'xticklabel', {},'yticklabel', {}, 'xlim', [-3, 9.5], 'ylim', [-1.5, 12.5] )
axis equal
xlabel( '$X$', 'fontsize', 40 )
ylabel( '$Y$', 'fontsize', 40 )

% ==================================== %
% Rhythmic DMP
% ==================================== %
tmp  = load( '../learned_parameters/rhythmic/heart.mat' );
data_rhythmic = tmp.data;

% The number of basis functions N
W_rhy = data_rhythmic.weight; 
g_rd  = data_rhythmic.goal;
[ ~, N ] = size( W_rhy );

% The gains for the transformation system
az = data_rhythmic.alpha_z; 
bz = data_rhythmic.beta_z;

% Duration of the discrete movement, demonstrated trajectory
tau_rd = data_rhythmic.tau;  

% The three elements of discrete DMP, subscript d is used
sr   = CanonicalSystem( 'rhythmic', tau_rd, 1.0 );
fr   = NonlinearForcingTerm( sr, N );
ts_r = TransformationSystem( az, bz, sr );

% The parameters for forward simulation
t0   = 0.0;
T     = (tau_rd*2*pi)*3;
dt    = 1e-4;
t_arr = 0:dt:T;
Nt    = length( t_arr );

% Calculate the nonlinear forcing term for rhythmic movement and rollout
input_arr = fr.calc_forcing_term( t_arr( 1:end-1 ), W_rhy, t0, eye( 2 ) );

% The input primitive, for rhythmic movement
prim_rhy  = input_arr + az * bz * g_rd;

[ y_arr2, ~, ~ ] = ts_r.rollout( zeros( 2, 1 ), zeros( 2, 1 ), zeros( 2, 1 ), prim_rhy, t0, t_arr  );

a2 = subplot( 3, 2, 2 );
plot( t_arr, sr.calc( t_arr ), 'linewidth', 4, 'color', c_orange);
set( a2,  'xlim', [0,3], 'xtick', 0:1:3, 'fontsize', 30, 'ylim', [ 0, 2*pi], 'ytick', 2*pi*[0., 0.5, 1.0], 'yticklabel', {'0', '$\pi$', '$2\pi$'} )
ylabel( '$s_r(t)$', 'fontsize', 40 )

a4 = subplot( 3, 2, 4 );
hold on
plot( a4, t_arr( 1:end-1 ), input_arr( 1, : ), 'linewidth', 4, 'linestyle', '-' , 'color', c_orange );
plot( a4, t_arr( 1:end-1 ), input_arr( 2, : ), 'linewidth', 4, 'linestyle', '-' , 'color', c_orange );
set( a4, 'yticklabel', {}, 'xlim', [0,3], 'xtick', 0:1:3, 'fontsize', 30 )
xlabel( '$t$', 'fontsize', 40 )
ylabel( '$\mathbf{f}_r(t)$', 'fontsize', 40 )

a6 = subplot( 3, 2, 6 );
hold on
plot( a6, y_arr2( 1, 1:end-1 ), y_arr2( 2, 1:end-1 ), 'linewidth', 4, 'color', c_orange );
scatter( a6, 0, 0, 100, 'filled',  'o', 'markerfacecolor', 'w', 'markeredgecolor', c_orange, 'linewidth', 4 )

set( a6, 'xticklabel', {}, 'yticklabel', {}, 'xlim', [-20,20], 'ylim',[-28,12])

xlabel( '$X$', 'fontsize', 40 )
ylabel( '$Y$', 'fontsize', 40 )
axis equal

fig_save( f, './images/fig1' )

%% (1Aa) [Figure 1a] Spatial Scaling

clear data*, close all; clc;

% ==================================== %
% Scaling for Discrete DMP
% ==================================== %
tmp  = load( '../learned_parameters/discrete/A_loose.mat' );
data_discrete = tmp.data;

% The number of basis functions N
W_dis = data_discrete.weight; 
g_dd  = data_discrete.goal;
[ ~, N ] = size( W_dis );

% The gains for the transformation system
as = data_discrete.alpha_s;
az = data_discrete.alpha_z; 
bz = data_discrete.beta_z;

% Duration of the discrete movement, demonstrated trajectory
tau_dd = data_discrete.tau;  

% The three elements of discrete DMP
sd   = CanonicalSystem( 'discrete', tau_dd, as );
fd   = NonlinearForcingTerm( sd, N );
ts_d = TransformationSystem( az, bz, sd );

f = figure( ); 
a1 = subplot( 1, 2, 1 );
hold on; axis equal
scl_arr    = [ 0.3, 0.5, 0.7, 0.9, 1.0, 1.2, 1.5, 2.0 ];
offset_arr = 10 * [ 0.2, 0.46, 0.9, 1.45, 2.2, 3.0, 4.0, 5.2];

% The Parameters for Forward Simulation
t0i   = 0.0;
T     = (data_discrete.tau)*2;
dt    = 1e-4;
t_arr = 0:dt:T;

for i = 1 : length( scl_arr )
    scl = scl_arr( i );
    off = offset_arr( i );

    % Calculate the nonlinear forcing term for discrete movement and rollout
    input_arr = fd.calc_forcing_term( t_arr( 1:end-1 ), W_dis, t0i, eye( 2 ) );
   
    % The input primitive, without scaling 
    prim_dis = input_arr + az * bz * g_dd;

    [ y_arr, ~, ~ ] = ts_d.rollout( zeros( 2, 1 ), zeros( 2, 1 ), zeros( 2, 1 ), scl*prim_dis, t0i, t_arr  );  

    % To emphasize
    if scl ~= 1
        plot( a1, off + y_arr( 1, : ), y_arr( 2,  : ), 'color', c_blue, 'linewidth', 3 )
    else
        plot( a1, off + y_arr( 1, : ), y_arr( 2,  : ), 'color', c_blue, 'linewidth', 8 )
    end

    scatter( a1, off + y_arr( 1,   1 ), y_arr( 2,   1 ), 300*scl, 'filled',  'o', 'markerfacecolor', 'w', 'markeredgecolor', c_blue, 'linewidth', 3 )
    scatter( a1, off + y_arr( 1, end ), y_arr( 2, end ), 300*scl, 'filled',  'd', 'markerfacecolor', 'w', 'markeredgecolor', c_blue, 'linewidth', 3 )

end
set( a1, 'xlim', [ 1, 67 ], 'ylim', [-25, 41], 'xticklabel', {}, 'yticklabel', {})
xlabel( 'X', 'fontsize', 40 )
ylabel( 'Y', 'fontsize', 40 )

% ==================================== %
% Scaling for Rhythmic DMP
% ==================================== %

a2 = subplot( 1, 2, 2 );
hold on; axis equal

% For rhythmic movement, load the data
tmp  = load( '../learned_parameters/rhythmic/heart.mat' );
data_rhythmic = tmp.data;

% The number of basis functions N
W_rhy = data_rhythmic.weight; 
g_rd  = data_rhythmic.goal;
[ ~, N ] = size( W_rhy );

% The gains for the transformation system
az = data_rhythmic.alpha_z; 
bz = data_rhythmic.beta_z;

% Duration of the discrete movement, demonstrated trajectory
tau_rd = data_rhythmic.tau;  

% The three elements of rhythmic DMP
% Note that this definition of the Canonical System is the simple linear function, 
% as we don't have a network of rhythmic Canonical System
sr   = CanonicalSystem( 'rhythmic', tau_rd, 1.0 );
fr   = NonlinearForcingTerm( sr, N );
ts_r = TransformationSystem( az, bz, sr );

% The Parameters for Forward Simulation
t0i   = 0.0;
T     = (tau_rd*2*pi)*2;
dt    = 1e-4;
t_arr = 0:dt:T;

for i = 1 : length( scl_arr )
    scl = scl_arr( i );

    % Calculate the nonlinear forcing term for discrete movement and rollout
    input_arr = fr.calc_forcing_term( t_arr( 1:end-1 ), W_rhy, t0i, eye( 2 ) );
   
    % Define input primitives 
    prim_rhy = input_arr + az*bz*g_rd;

    % Rollout
    [ y_arr, ~, ~ ] = ts_r.rollout( zeros( 2, 1 ), zeros( 2, 1 ), zeros( 2, 1 ), scl*prim_rhy, t0i, t_arr  );  

    if scl ~= 1
        plot( a2, y_arr( 1, : ), y_arr( 2,  : ), 'color', c_orange, 'linewidth', 3 )
    else
        plot( a2, y_arr( 1, : ), y_arr( 2,  : ), 'color', c_orange, 'linewidth', 8 )
    end
end

scatter( a2, 0, 0, 100, 'filled',  'o', 'markerfacecolor', 'w', 'markeredgecolor', c_orange, 'linewidth', 3 )
xlabel( 'X', 'fontsize', 40 )
ylabel( 'Y', 'fontsize', 40 )

set( a2, 'xlim', [-36, 36 ], 'ylim', [-51, 21], 'xticklabel', {}, 'yticklabel', {} )
fig_save( f, './images/fig1a' )


%% (1Ab) [Figure 1b] Temporal Scaling 

clear data*, close all; clc;

% ==================================== %
% Temporal Scaling for Discrete DMP
% ==================================== %
f = figure( );

% Load the A alphabet 
tmp  = load( '../learned_parameters/discrete/A_loose.mat' );
data_discrete = tmp.data;

% The number of basis functions N
W_dis = data_discrete.weight; 
g_dd  = data_discrete.goal;
[ ~, N ] = size( W_dis );

% The gains for the transformation system
as = data_discrete.alpha_s;
az = data_discrete.alpha_z; 
bz = data_discrete.beta_z;

% Duration of the discrete movement, demonstrated trajectory
tau_dd = data_discrete.tau;  

% The three elements of discrete DMP
sd   = CanonicalSystem( 'discrete', tau_dd, as );
fd   = NonlinearForcingTerm( sd, N );
ts_d = TransformationSystem( az, bz, sd );

% First, generate the trajectories
ts_arr = [ 0.5, 1.0, 20 ];

% The Parameters for Forward Simulation
t0i   = 0.0;
T     = tau_dd*2;
dt    = 1e-3;
t_arr = 0:dt:(T*20);

% The generated trajectory
y_data_arr = zeros( 2, length( t_arr ), length( ts_arr ) );

for i = 1 : length( ts_arr )
    ts = ts_arr( i );

    % Calculate the nonlinear forcing term for discrete movement and rollout
    input_arr = fd.calc_forcing_term( t_arr( 1:end-1 )/ts, W_dis, t0i, eye( 2 ) );

    % The input primitives
    prim_arr = input_arr + az * bz * g_dd;
  
    [ y_arr, ~, ~ ] = ts_d.rollout( zeros( 2, 1 ), zeros( 2, 1 ), zeros( 2, 1), prim_arr, t0i, t_arr  );  

    y_data_arr( :, :, i ) = y_arr;
    
end

ttmp_arr = tau_dd *(0.25:0.25:2.00);
offset   = 10*( 1:length( ttmp_arr ) );

% Going through the iterations
% Over the time array
for i = 1: length( ts_arr )

    a = subplot( length( ts_arr ), 2, 2*i-1, 'parent', f );
    hold on; axis equal

    % Get the trajectory data
    y_tmp = y_data_arr( :, :, i );

    for j = 1 : length( ttmp_arr )
        tt = ttmp_arr( j );

        [ ~, idx] = min( abs( t_arr - tt ) );
        off = offset( j );

        plot( a, off + y_tmp( 1, 1:idx ), y_tmp( 2, 1:idx ), 'linewidth', 5, 'color', c_blue )
        scatter( a, off + y_tmp( 1, idx ), y_tmp( 2, idx ), 100, 'markerfacecolor', 'w', 'markeredgecolor', c_blue, 'linewidth', 3 )        
    end
    set( a, 'xlim', [3, max(offset)+10], 'xticklabel', {}, 'yticklabel', {}  )
    
    title( a, ['$s_t=', num2str( ts_arr( i ), '%.2f' ), '$' ], 'fontsize', 40 )
end

% ==================================== %
% Temporal Scaling for Rhythmic DMP
% ==================================== %

% For rhythmic movement, load the data
tmp  = load( '../learned_parameters/rhythmic/heart.mat' );
data_rhythmic = tmp.data;

% The number of basis functions N
W_rhy = data_rhythmic.weight; 
g_rd  = data_rhythmic.goal;
[ ~, N ] = size( W_rhy );

% The gains for the transformation system
az = data_rhythmic.alpha_z; 
bz = data_rhythmic.beta_z;

% Duration of the discrete movement, demonstrated trajectory
tau_rd = data_rhythmic.tau;  

% The three elements of rhythmic DMP
% Note that this definition of the Canonical System is the simple linear function, 
% as we don't have a network of rhythmic Canonical System
sr   = CanonicalSystem( 'rhythmic', tau_rd, 1.0 );
fr   = NonlinearForcingTerm( sr, N );
ts_r = TransformationSystem( az, bz, sr );

% First, generate the trajectories
tscl_arr = [ 0.5, 1.0, 2.0 ];

% The Parameters for Forward Simulation
t0i   = 0.0;
T     = (tau_rd*2*pi)*2;
dt    = 1e-3;
t_arr = 0:dt:T;

y_data_arr = zeros( 2, length( t_arr ), length( tscl_arr ) );
scl = 1.8;

for i = 1 : length( ts_arr )
    ts = ts_arr( i );

    % Calculate the nonlinear forcing term for discrete movement and rollout
    input_arr = fr.calc_forcing_term( t_arr( 1:end-1 )/ts, W_rhy, t0i, eye( 2 ) );

    % The input primitives
    prim_rhy = input_arr + az * bz * g_rd;

    [ y_arr, ~, ~ ] = ts_r.rollout( zeros( 2, 1), zeros( 2, 1 ), zeros( 2, 1 ), scl*prim_rhy, t0i, t_arr  );  

    y_data_arr( :, :, i ) = y_arr;
    
end

ttmp_arr = tau_rd * (0.25:0.25:2.0)*2*pi;
offset  = 66*( 1:length( ttmp_arr ) );

% Going through the iterations
% Over the time array
for i = 1: length( ts_arr )

    a = subplot( length( ts_arr ), 2, 2*i, 'parent', f );
    hold on; axis equal
    y_tmp = y_data_arr( :, :, i );
    
    for j = 1 : length( ttmp_arr )
        tt = ttmp_arr( j );

        [ ~, idx] = min( abs( t_arr - tt ) );
        off = offset( j );
        plot( a, off + y_tmp( 1, 1:idx ), y_tmp( 2, 1:idx ), 'linewidth', 5, 'color', c_orange )
        scatter( a, off + y_tmp( 1, idx ), y_tmp( 2, idx ), 100, 'markerfacecolor', 'w', 'markeredgecolor', c_orange, 'linewidth', 3 )
    end

    set( a, 'xlim', [30, max(offset)+45], 'xticklabel', {}, 'yticklabel', {}  )
    title( a, ['$s_t=', num2str( tscl_arr( i ), '%.2f' ), '$'], 'fontsize', 40 )
    
end

fig_save( f, './images/fig1b' )

%% (1Ac) [Figure 1c] Rotational Scaling

clear data*, close all; clc;

% ==================================== %
% Rotation for Discrete DMP
% ==================================== %
% Load the A alphabet 
tmp  = load( '../learned_parameters/discrete/A_loose.mat' );
data_discrete = tmp.data;

% The number of basis functions N
W_dis = data_discrete.weight; 
g_dd  = data_discrete.goal;
[ ~, N ] = size( W_dis );

% The gains for the transformation system
as = data_discrete.alpha_s;
az = data_discrete.alpha_z; 
bz = data_discrete.beta_z;

% Duration of the discrete movement, demonstrated trajectory
tau_dd = data_discrete.tau;  

% The three elements of discrete DMP
sd   = CanonicalSystem( 'discrete', tau_dd, as );
fd   = NonlinearForcingTerm( sd, N );
ts_d = TransformationSystem( az, bz, sd );

f = figure( ); 
a1 = subplot( 1, 2, 1 );
hold on; axis equal

% The angular displacement
ang_arr = linspace( 0, 2*pi, 6 );
ang_arr = ang_arr( 1:end-1 );

% The Parameters for Forward Simulation
t0i   = 0.0;
T     = tau_dd*2;
dt    = 1e-4;
t_arr = 0:dt:T;
scl   = 2.5;

for i = 1 : length( ang_arr )

    ang = ang_arr( i );

    % Calculate the nonlinear forcing term for discrete movement and rollout
    input_arr = fd.calc_forcing_term( t_arr( 1:end-1 ), W_dis, t0i, eye( 2 ) );
 
    % Rotation matrices 
    R = [ cos( ang ), -sin( ang );
          sin( ang ),  cos( ang )];

    % The input primitive
    prim_dis = scl*R*( input_arr + az*bz*g_dd );

    [ y_arr, ~, ~ ] = ts_d.rollout( zeros( 2, 1 ), zeros( 2, 1 ), zeros( 2, 1 ), prim_dis, t0i, t_arr  );  
    
    % The offset 
    y0 = R * [ 10; 0 ];

    plot( a1, y0( 1 ) + y_arr( 1, : ), y0( 2 ) + y_arr( 2,  : ), 'color', c_blue )
    scatter( a1, y0( 1 ) +y_arr( 1,   1 ),y0( 2 ) + y_arr( 2,   1 ), 200, 'filled',  'o', 'markerfacecolor', 'w', 'markeredgecolor', c_blue, 'linewidth', 3 )
    scatter( a1, y0( 1 ) +y_arr( 1, end ), y0( 2 ) +y_arr( 2, end ), 200, 'filled',  'd', 'markerfacecolor', 'w', 'markeredgecolor', c_blue, 'linewidth', 3 )

end
set( a1, 'xlim', get( a1, 'xlim' ) + [-4, 2], 'ylim', get( a1, 'ylim' ) + [-3, 3], 'xticklabel',{},'yticklabel',{} )
xlabel( a1, 'X', 'fontsize', 40)
ylabel( a1, 'Y', 'fontsize', 40)

% ==================================== %
% Rotation for Rhythmic DMP
% ==================================== %
a2 = subplot( 1, 2, 2 );
hold on; axis equal

% For rhythmic movement, load the data
tmp  = load( '../learned_parameters/rhythmic/heart.mat' );
data_rhythmic = tmp.data;

% The number of basis functions N
W_rhy = data_rhythmic.weight; 
g_rd  = data_rhythmic.goal;
[ ~, N ] = size( W_rhy );

% The gains for the transformation system
az = data_rhythmic.alpha_z; 
bz = data_rhythmic.beta_z;

% Duration of the discrete movement, demonstrated trajectory
tau_rd = data_rhythmic.tau;  

% The three elements of rhythmic DMP
% Note that this definition of the Canonical System is the simple linear function, 
% as we don't have a network of rhythmic Canonical System
sr   = CanonicalSystem( 'rhythmic', tau_rd, 1.0 );
fr   = NonlinearForcingTerm( sr, N );
ts_r = TransformationSystem( az, bz, sr );

% The Parameters for Forward Simulation
t0i   = 0.0;
T     = (tau_rd*2*pi)*2;
dt    = 1e-4;
t_arr = 0:dt:T;

scl = 0.6;

for i = 1 : length( ang_arr )

    ang = ang_arr( i );

    % Calculate the nonlinear forcing term for discrete movement and rollout
    input_arr = fr.calc_forcing_term( t_arr( 1:end-1 ), W_rhy, t0i, eye( 2 ) );
 
    % Rotation matrices 
    R = [ cos( ang ), -sin( ang );
          sin( ang ),  cos( ang )];

    % The input primitive
    prim_rhy = scl*R*( input_arr + az*bz*g_rd );

    [ y_arr, ~, ~ ] = ts_r.rollout( zeros( 2, 1 ), zeros( 2, 1 ), zeros( 2, 1 ), prim_rhy, t0i, t_arr  );  

    % The offset
    y0 = R*[15;0];    

    plot( a2, y0( 1 ) + y_arr( 1, : ), y0( 2 ) + y_arr( 2,  : ), 'color', c_orange )
    scatter( a2, y0( 1 ),y0( 2 ), 200, 'filled',  'o', 'markerfacecolor', 'w', 'markeredgecolor', c_orange, 'linewidth', 3 )
    
end

set( a2, 'xlim', get( a2, 'xlim' ) + [-4, 2], 'ylim', get( a2, 'ylim' ) + [-3, 3], 'xticklabel',{},'yticklabel',{} )
xlabel( a2, 'X', 'fontsize', 40)
ylabel( a2, 'Y', 'fontsize', 40)

fig_save( f, './images/fig1c' )

%% (1B) [Figure 2a] Merging Two Discrete Movements

close all; clc;

f = figure( ); 

% ==================================== %
% Combination of two Discrete Movements
% For simple combination
% ==================================== %

% Load two mat functions, discrete movements
tmp1 = load( '../learned_parameters/discrete/A_loose.mat' ); dis1 = tmp1.data;
tmp2 = load( '../learned_parameters/discrete/B_loose.mat' ); dis2 = tmp2.data;

% Get the degrees of freedom and the number of weights 
n1 = size( dis1.weight, 1 ); N1 = size( dis1.weight, 2 );
n2 = size( dis2.weight, 1 ); N2 = size( dis2.weight, 2 );

assert( n1 == n2 && N1 == N2 )
n = n1; N = N1;

% The alpha, beta values must be same
% The gains for the transformation system
as1 = dis1.alpha_s;   as2 = dis2.alpha_s;
az1 = dis1.alpha_z;   az2 = dis2.alpha_z; 
bz1 = dis1.beta_z;    bz2 = dis2.beta_z;

assert( as1 == as2 && az1 == az2 && bz1 == bz2 );
az = az1; as = as1; bz = bz1;

g_dd1   = dis1.goal; g_dd2   = dis2.goal;
tau_dd1 = dis1.tau;  tau_dd2 = dis2.tau;

tau = tau_dd1;

% The DMP for the input
sd  = CanonicalSystem( 'discrete', tau, as );
ts  = TransformationSystem( az, bz, sd );

% Rollout for each discrete movement
t0i   =    0.0;    % The initial time of the movement rollout
dt    =   1e-4;    % Time-step  for the simulation
t_arr = 0:dt:(tau*6);  % Time array for the simulation

% The discrete movements' DMPs
sd_1  = CanonicalSystem( 'discrete', tau_dd1, as );
sd_2  = CanonicalSystem( 'discrete', tau_dd2, as );

fd_1  = NonlinearForcingTerm( sd_1, N );
fd_2  = NonlinearForcingTerm( sd_2, N );

% Time scaling factor
tmp1 = 1.0; tmp2 = 1.0;

tscl1 = tau/tau_dd1 * tmp1;
tscl2 = tau/tau_dd2 * tmp2;

input_arr1 = fd_1.calc_forcing_term( t_arr( 1:end-1 )/tscl1, dis1.weight, t0i, eye( 2 ) );
input_arr2 = fd_2.calc_forcing_term( t_arr( 1:end-1 )/tscl2, dis2.weight, t0i, eye( 2 ) );

prim1 = input_arr1 + az*bz*g_dd1;
prim2 = input_arr2 + az*bz*g_dd2;

a = subplot( 3, 1, 1, 'parent', f );
hold on; axis equal

gain_arr = 0:0.1:1;
Ntmp = length( gain_arr );
off = 12.0;

scl1 = 1.1;
scl2 = 0.9;

tscl_arr = linspace( 1, 3, Ntmp );

for i = 1: Ntmp
    gain = gain_arr( Ntmp-i+1 );

    prim = gain * prim1 * scl1 + ( 1 - gain ) * prim2 * scl2;

    [ y_arr_comb, ~, ~] = ts.rollout( zeros( n, 1 ), zeros( n, 1 ), zeros( n, 1 ), prim, t0i, t_arr  );    
    
    plot( a, (i-1)*off+y_arr_comb( 1, : )-3, (i-1)*0.08*off+y_arr_comb( 2, : )+3.5, 'linewidth', 6, 'color', c_blue )
end
set( a, 'xticklabel', {}, 'yticklabel', {}, 'xlim', [-12, 132], 'ylim', [0, 20]  )
set( a, 'color', 'none' )

% ==================================== %
% Combination of two Discrete Movements
% For Temporal Scaling
% ==================================== %

a = subplot( 3, 1, 2, 'parent', f );
hold on; axis equal

gain_arr = 0:0.1:1;
Ntmp = length( gain_arr );
off = 12.0;

scl1 = 1.1;
scl2 = 0.9;

% Time scaling factor
tmp1 = 2.0; tmp2 = 1.0;

tscl1 = tau/tau_dd1 * tmp1;
tscl2 = tau/tau_dd2 * tmp2;

input_arr1 = fd_1.calc_forcing_term( t_arr( 1:end-1 )/tscl1, dis1.weight, t0i, eye( 2 ) );
input_arr2 = fd_2.calc_forcing_term( t_arr( 1:end-1 )/tscl2, dis2.weight, t0i, eye( 2 ) );

prim1 = input_arr1 + az*bz*g_dd1;
prim2 = input_arr2 + az*bz*g_dd2;


for i = 1: Ntmp
    gain = gain_arr( Ntmp-i+1 );
    prim = gain * prim1 * scl1 + ( 1 - gain ) * prim2 * scl2;

    [ y_arr_comb, ~, ~] = ts.rollout( zeros( n, 1 ), zeros( n, 1 ), zeros( n, 1 ), prim, t0i, t_arr  );    
    
    plot( a, (i-1)*off+y_arr_comb( 1, : )-3, (i-1)*0.08*off+y_arr_comb( 2, : )+3.5, 'linewidth', 6, 'color', c_blue )
end
set( a, 'xticklabel', {}, 'yticklabel', {}, 'xlim', [-12, 132], 'ylim', [0, 20]  )
set( a, 'color', 'none' )


% ==================================== %
% Combination of two Discrete Movements
% For Spatial Scaling
% ==================================== %

a = subplot( 3, 1, 3, 'parent', f );
hold on; axis equal

gain_arr = 0:0.1:1;
Ntmp = length( gain_arr );
off = 12.0;

scl1 = 1.1;
scl2 = 0.4;

% Time scaling factor
tmp1 = 1.0; tmp2 = 1.0;

tscl1 = tau/tau_dd1 * tmp1;
tscl2 = tau/tau_dd2 * tmp2;

input_arr1 = fd_1.calc_forcing_term( t_arr( 1:end-1 )/tscl1, dis1.weight, t0i, eye( 2 ) );
input_arr2 = fd_2.calc_forcing_term( t_arr( 1:end-1 )/tscl2, dis2.weight, t0i, eye( 2 ) );

prim1 = input_arr1 + az*bz*g_dd1;
prim2 = input_arr2 + az*bz*g_dd2;

S = [ cos( pi/4 ), -sin( pi/4 ); 
      sin( pi/4 ),  cos( pi/4 ) ];

for i = 1: Ntmp
    gain = gain_arr( Ntmp-i+1 );
    prim = S * gain * prim1 * scl1 + inv( S ) * ( 1 - gain ) * prim2 * scl2;

    [ y_arr_comb, ~, ~] = ts.rollout( zeros( n, 1 ), zeros( n, 1 ), zeros( n, 1 ), prim, t0i, t_arr  );    
    
    plot( a, (i-1)*off+y_arr_comb( 1, : ), (i-1)*0.09*off+y_arr_comb( 2, : )+3.5, 'linewidth', 6, 'color', c_blue )
end
set( a, 'xticklabel', {}, 'yticklabel', {}, 'xlim', [-12, 132], 'ylim', [0, 20]  )
set( a, 'color', 'none' )

fig_save( f, './images/fig2a_D' )

%% (1B) [Figure 2b] Merging Discrete and Rhythmic Movements

close all; clc;

f = figure( ); 

% ==================================== %
% Combination of Discrete and Rhytmic Movements
% For simple combination
% ==================================== %

% Load two mat functions, discrete movements
tmp1 = load( '../learned_parameters/discrete/A_loose.mat' ); dis1 = tmp1.data;
tmp2 = load( '../learned_parameters/rhythmic/heart.mat' );   rhy1 = tmp2.data;

% Get the degrees of freedom and the number of weights 
n1 = size( dis1.weight, 1 ); N1 = size( dis1.weight, 2 );
n2 = size( rhy1.weight, 1 ); N2 = size( rhy1.weight, 2 );

assert( n1 == n2 && N1 == N2 )
n = n1; N = N1;

% The alpha, beta values must be same
% The gains for the transformation system
as1 = dis1.alpha_s;   as2 = dis1.alpha_s;
az1 = dis1.alpha_z;   az2 = rhy1.alpha_z; 
bz1 = dis1.beta_z;    bz2 = rhy1.beta_z;

assert( as1 == as2 && az1 == az2 && bz1 == bz2 );
az = az1; as = as1; bz = bz1;

g_dd   = dis1.goal; g_rd   = rhy1.goal;
tau_dd = dis1.tau;  tau_rd = rhy1.tau;

tau = tau_rd;

% The DMP for the input
% Either discrete or rhythmic can be used.
s  = CanonicalSystem( 'rhythmic', tau, as );
ts = TransformationSystem( az, bz, s );

% Rollout for each discrete movement
t0i   =    0.0;    % The initial time of the movement rollout
dt    =   1e-4;    % Time-step  for the simulation
t_arr = 0:dt:(tau_rd*2*pi*3);  % Time array for the simulation

% The discrete movements' DMPs
sd  = CanonicalSystem( 'discrete', tau_dd, as );
sr  = CanonicalSystem( 'rhythmic', tau_rd, as );

fd  = NonlinearForcingTerm( sd, N );
fr  = NonlinearForcingTerm( sr, N );

% Time scaling factor
tmp1 = 10.0; tmp2 = 1.0;

tscl1 = tau/tau_dd * tmp1;
tscl2 = tau/tau_rd * tmp2;

input_arr1 = fd.calc_forcing_term( t_arr( 1:end-1 )/tscl1, dis1.weight, t0i, eye( 2 ) );
input_arr2 = fr.calc_forcing_term( t_arr( 1:end-1 )/tscl2, rhy1.weight, t0i, eye( 2 ) );

prim1 = input_arr1 + az*bz*g_dd;
prim2 = input_arr2 + az*bz*g_rd;

a = subplot( 3, 1, 1, 'parent', f );
hold on; axis equal

gain_arr = 0:0.1:1;
Ntmp = length( gain_arr );
off = 13.5;

scl1 =  1.1;
scl2 =  0.4;
tscl_arr = linspace( 1, 3, Ntmp );

for i = 1: Ntmp
    gain = gain_arr( Ntmp-i+1 );

    prim = gain * prim1 * scl1 + ( 1 - gain ) * prim2 * scl2;

    [ y_arr_comb, ~, ~] = ts.rollout( zeros( n, 1 ), zeros( n, 1 ), zeros( n, 1 ), prim, t0i, t_arr  );    
    if i == Ntmp
        plot( a, (i-1)*off+2+y_arr_comb( 1, : )-3, (i-1)*0.075*off+y_arr_comb( 2, : )+3.5, 'linewidth', 6, 'color', c_blue * gain + c_orange * (1-gain) )
    else
        plot( a, (i-1)*off+y_arr_comb( 1, : )-3, (i-1)*0.075*off+y_arr_comb( 2, : )+3.5, 'linewidth', 6, 'color', c_blue * gain + c_orange * (1-gain) )
    end
end
set( a, 'xticklabel', {}, 'yticklabel', {}, 'xlim', [-12, 147], 'ylim', [0, 20]  )
set( a, 'color', 'none' )


% ==================================== %
% Combination of Discrete and Rhythmic Movements
% For Temporal Scaling
% ==================================== %

a = subplot( 3, 1, 2, 'parent', f );
hold on; axis equal         

% Rollout for each discrete movement
t0i   =    0.0;    % The initial time of the movement rollout
dt    =   1e-4;    % Time-step  for the simulation
t_arr = 0:dt:(tau_rd*2*pi*40);  % Time array for the simulation

% Time scaling factor
tmp1 = 40.0; tmp2 = 1.0;

tscl1 = tau/tau_dd * tmp1;
tscl2 = tau/tau_rd * tmp2;

input_arr1 = fd.calc_forcing_term( t_arr( 1:end-1 )/tscl1, dis1.weight, t0i, eye( 2 ) );
input_arr2 = fr.calc_forcing_term( t_arr( 1:end-1 )/tscl2, rhy1.weight, t0i, eye( 2 ) );

prim1 = input_arr1 + az*bz*g_dd;
prim2 = input_arr2 + az*bz*g_rd;

a = subplot( 3, 1, 2, 'parent', f );
hold on; axis equal

gain_arr = 0:0.1:1;
Ntmp = length( gain_arr );
off = 13.5;

scl1 =  1.1;
scl2 =  0.4;
tscl_arr = linspace( 1, 3, Ntmp );

for i = 1: Ntmp
    gain = gain_arr( Ntmp-i+1 );
    prim = gain * prim1 * scl1 + ( 1 - gain ) * prim2 * scl2;

    [ y_arr_comb, ~, ~] = ts.rollout( zeros( n, 1 ), zeros( n, 1 ), zeros( n, 1 ), prim, t0i, t_arr  );    
    
    if i == Ntmp
        plot( a, (i-1)*off+2+y_arr_comb( 1, : )-3, (i-1)*0.075*off+y_arr_comb( 2, : )+3.5, 'linewidth', 6, 'color', c_blue * gain + c_orange * (1-gain) )
    else
        plot( a, (i-1)*off+y_arr_comb( 1, : )-3, (i-1)*0.075*off+y_arr_comb( 2, : )+3.5, 'linewidth', 6, 'color', c_blue * gain + c_orange * (1-gain) )
    end

end
set( a, 'xticklabel', {}, 'yticklabel', {}, 'xlim', [-12, 145], 'ylim', [0, 20]  )
set( a, 'color', 'none' )

% ==================================== %
% Combination of Discrete and Rhythmic Movements
% For Spatial Scaling
% ==================================== %

a = subplot( 3, 1, 3, 'parent', f );
hold on; axis equal         

% Rollout for each discrete movement
t0i   =    0.0;    % The initial time of the movement rollout
dt    =   1e-4;    % Time-step  for the simulation
t_arr = 0:dt:(tau_rd*2*pi*40);  % Time array for the simulation

% Time scaling factor
tmp1 = 40.0; tmp2 = 1.0;

tscl1 = tau/tau_dd * tmp1;
tscl2 = tau/tau_rd * tmp2;

input_arr1 = fd.calc_forcing_term( t_arr( 1:end-1 )/tscl1, dis1.weight, t0i, eye( 2 ) );
input_arr2 = fr.calc_forcing_term( t_arr( 1:end-1 )/tscl2, rhy1.weight, t0i, eye( 2 ) );

prim1 = input_arr1 + az*bz*g_dd;
prim2 = input_arr2 + az*bz*g_rd;

a = subplot( 3, 1, 3, 'parent', f );
hold on; axis equal

gain_arr = 0:0.1:1;
Ntmp = length( gain_arr );
off = 13.5;

scl1 =  1.1;
scl2 =  0.2;
tscl_arr = linspace( 1, 3, Ntmp );

S = [ cos( pi/4 ), -sin( pi/4 ); 
      sin( pi/4 ),  cos( pi/4 ) ];

for i = 1: Ntmp
    gain = gain_arr( Ntmp-i+1 );
    prim = S*gain * prim1 * scl1 + inv( S ) * ( 1 - gain ) * prim2 * scl2;

    [ y_arr_comb, ~, ~] = ts.rollout( zeros( n, 1 ), zeros( n, 1 ), zeros( n, 1 ), prim, t0i, t_arr  );    
    
    if i == Ntmp
        plot( a, (i-1)*off+2+y_arr_comb( 1, : )-3, (i-1)*0.070*off+y_arr_comb( 2, : )+3.5, 'linewidth', 6, 'color', c_blue * gain + c_orange * (1-gain) )
    else
        plot( a, (i-1)*off+y_arr_comb( 1, : )-3, (i-1)*0.070*off+y_arr_comb( 2, : )+3.5, 'linewidth', 6, 'color', c_blue * gain + c_orange * (1-gain) )
    end
    
end
set( a, 'xticklabel', {}, 'yticklabel', {}, 'xlim', [-12, 145], 'ylim', [0, 20]  )
set( a, 'color', 'none' )

fig_save( f, './images/fig2a_E' )

%% (1B) [Figure 2c] Merging Two Rhythmic Movements

close all; clc;

f = figure( ); 

% ==================================== %
% Combination of Discrete and Rhytmic Movements
% For simple combination
% ==================================== %

% Load two mat functions, discrete movements
tmp1 = load( '../learned_parameters/rhythmic/heart.mat' );  rhy1 = tmp1.data;
tmp2 = load( '../learned_parameters/rhythmic/circle.mat' ); rhy2 = tmp2.data;

% Get the degrees of freedom and the number of weights 
n1 = size( rhy1.weight, 1 ); N1 = size( rhy1.weight, 2 );
n2 = size( rhy2.weight, 1 ); N2 = size( rhy2.weight, 2 );

assert( n1 == n2 && N1 == N2 )
n = n1; N = N1;

% The alpha, beta values must be same
% The gains for the transformation system
az1 = rhy1.alpha_z;   az2 = rhy2.alpha_z; 
bz1 = rhy1.beta_z;    bz2 = rhy2.beta_z;

assert( az1 == az2 && bz1 == bz2 );
az = az1; as = as1; bz = bz1;

g_rd1   = rhy1.goal; g_rd2   = rhy2.goal;
tau_rd1 = rhy1.tau;  tau_rd2 = rhy2.tau;

tau = tau_rd1;

% The DMP for the input
% Either discrete or rhythmic can be used.
s  = CanonicalSystem( 'rhythmic', tau, as );
ts = TransformationSystem( az, bz, s );

% Rollout for each discrete movement
t0i   =    0.0;    % The initial time of the movement rollout
dt    =   1e-4;    % Time-step  for the simulation
t_arr = 0:dt:(tau_rd1*2*pi*3);  % Time array for the simulation

% The discrete movements' DMPs
sr1  = CanonicalSystem( 'rhythmic', tau_rd1, 1.0 );
sr2  = CanonicalSystem( 'rhythmic', tau_rd2, 1.0 );

fr1  = NonlinearForcingTerm( sr1, N );
fr2  = NonlinearForcingTerm( sr2, N );

% Time scaling factor
tmp1 = 1.0; tmp2 = 1.0;

tscl1 = tau/tau_rd1 * tmp1;
tscl2 = tau/tau_rd2 * tmp2;

input_arr1 = fr1.calc_forcing_term( t_arr( 1:end-1 )/tscl1, rhy1.weight, t0i, eye( 2 ) );
input_arr2 = fr2.calc_forcing_term( t_arr( 1:end-1 )/tscl2, rhy2.weight, t0i, eye( 2 ) );

prim1 = input_arr1 + az*bz*g_rd1;
prim2 = input_arr2 + az*bz*g_rd2;

a = subplot( 3, 1, 1, 'parent', f );
hold on; axis equal

gain_arr = 0:0.1:1;
Ntmp = length( gain_arr );
off = 14.0;

scl1 =  0.4;
scl2 = 11.0;
tscl_arr = linspace( 1, 3, Ntmp );

for i = 1: Ntmp
    gain = gain_arr( Ntmp-i+1 );

    prim = gain * prim1 * scl1 + ( 1 - gain ) * prim2 * scl2;

    [ y_arr_comb, ~, ~] = ts.rollout( zeros( n, 1 ), zeros( n, 1 ), zeros( n, 1 ), prim, t0i, t_arr  );    

    plot( a, (i-1)*off+y_arr_comb( 1, : )-3, -(i-1)*0.020*off+y_arr_comb( 2, : )+3.5, 'linewidth', 3, 'color', c_orange )

end
set( a, 'xticklabel', {}, 'yticklabel', {}, 'xlim', [-15, 142], 'ylim', [-10, 10]  )
set( a, 'color', 'none' )


% Rollout for each discrete movement
t0i   =    0.0;    % The initial time of the movement rollout
dt    =   1e-4;    % Time-step  for the simulation
t_arr = 0:dt:(tau_rd1*2*pi*20);  % Time array for the simulation


% Time scaling factor
tmp1 = 1.0; tmp2 = 10.0;

tscl1 = tau/tau_rd1 * tmp1;
tscl2 = tau/tau_rd2 * tmp2;

input_arr1 = fr1.calc_forcing_term( t_arr( 1:end-1 )/tscl1, rhy1.weight, t0i, eye( 2 ) );
input_arr2 = fr2.calc_forcing_term( t_arr( 1:end-1 )/tscl2, rhy2.weight, t0i, eye( 2 ) );

prim1 = input_arr1 + az*bz*g_rd1;
prim2 = input_arr2 + az*bz*g_rd2;

a = subplot( 3, 1, 2, 'parent', f );
hold on; axis equal

gain_arr = 0:0.1:1;
Ntmp = length( gain_arr );
off = 14.0;

scl1 =  0.4;
scl2 = 11.0;
tscl_arr = linspace( 1, 3, Ntmp );

for i = 1: Ntmp
    gain = gain_arr( Ntmp-i+1 );

    prim = gain * prim1 * scl1 + ( 1 - gain ) * prim2 * scl2;

    [ y_arr_comb, ~, ~] = ts.rollout( zeros( n, 1 ), zeros( n, 1 ), zeros( n, 1 ), prim, t0i, t_arr  );    

    plot( a, (i-1)*off+y_arr_comb( 1, : )-3, -(i-1)*0.020*off+y_arr_comb( 2, : )+3.5, 'linewidth', 3, 'color', c_orange )

end
set( a, 'xticklabel', {}, 'yticklabel', {}, 'xlim', [-15, 142], 'ylim', [-10, 10]  )
set( a, 'color', 'none' )



% Rollout for each discrete movement
t0i   =    0.0;    % The initial time of the movement rollout
dt    =   1e-4;    % Time-step  for the simulation
t_arr = 0:dt:(tau_rd1*2*pi*20);  % Time array for the simulation


% Time scaling factor
tmp1 = 1.0; tmp2 = 10.0;

tscl1 = tau/tau_rd1 * tmp1;
tscl2 = tau/tau_rd2 * tmp2;

input_arr1 = fr1.calc_forcing_term( t_arr( 1:end-1 )/tscl1, rhy1.weight, t0i, eye( 2 ) );
input_arr2 = fr2.calc_forcing_term( t_arr( 1:end-1 )/tscl2, rhy2.weight, t0i, eye( 2 ) );

prim1 = input_arr1 + az*bz*g_rd1;
prim2 = input_arr2 + az*bz*g_rd2;

a = subplot( 3, 1, 3, 'parent', f );
hold on; axis equal

gain_arr = 0:0.1:1;
Ntmp = length( gain_arr );
off = 14.0;

scl1 =  0.4;
scl2 = 11.0;
tscl_arr = linspace( 1, 3, Ntmp );


S = [ cos( pi/3 ), -sin( pi/3 ); 
      sin( pi/3 ),  cos( pi/3 ) ];


for i = 1: Ntmp
    gain = gain_arr( Ntmp-i+1 );

    prim = S^2 * gain * prim1 * scl1 + ( 1 - gain ) * prim2 * scl2;

    [ y_arr_comb, ~, ~] = ts.rollout( zeros( n, 1 ), zeros( n, 1 ), zeros( n, 1 ), prim, t0i, t_arr  );    

    plot( a, (i-1)*off+y_arr_comb( 1, : )-3, (i-1)*0.020*off+y_arr_comb( 2, : )+3.5, 'linewidth', 3, 'color', c_orange )

end
set( a, 'xticklabel', {}, 'yticklabel', {}, 'xlim', [-15, 142], 'ylim', [-5, 15]  )
set( a, 'color', 'none' )



fig_save( f, './images/fig2a_F' )

%% (1B) [Figure 2c] Merging Two Rhythmic Movements

close all; clc;

f = figure( ); 

% ==================================== %
% Combination of Discrete and Rhytmic Movements
% For simple combination
% ==================================== %

% Load two mat functions, discrete movements
tmp1 = load( '../learned_parameters/rhythmic/heart.mat' );  rhy1 = tmp1.data;
tmp2 = load( '../learned_parameters/rhythmic/circle.mat' ); rhy2 = tmp2.data;

% Get the degrees of freedom and the number of weights 
n1 = size( rhy1.weight, 1 ); N1 = size( rhy1.weight, 2 );
n2 = size( rhy2.weight, 1 ); N2 = size( rhy2.weight, 2 );

assert( n1 == n2 && N1 == N2 )
n = n1; N = N1;

% The alpha, beta values must be same
% The gains for the transformation system
az1 = rhy1.alpha_z;   az2 = rhy2.alpha_z; 
bz1 = rhy1.beta_z;    bz2 = rhy2.beta_z;

assert( az1 == az2 && bz1 == bz2 );
az = az1; as = as1; bz = bz1;

g_rd1   = rhy1.goal; g_rd2   = rhy2.goal;
tau_rd1 = rhy1.tau;  tau_rd2 = rhy2.tau;

tau = tau_rd1;

% The DMP for the input
% Either discrete or rhythmic can be used.
s  = CanonicalSystem( 'rhythmic', tau, as );
ts = TransformationSystem( az, bz, s );

% Rollout for each discrete movement
t0i   =    0.0;    % The initial time of the movement rollout
dt    =   1e-4;    % Time-step  for the simulation
t_arr = 0:dt:(tau_rd1*2*pi*10);  % Time array for the simulation

% The discrete movements' DMPs
sr1  = CanonicalSystem( 'rhythmic', tau_rd1, 1.0 );
sr2  = CanonicalSystem( 'rhythmic', tau_rd2, 1.0 );

fr1  = NonlinearForcingTerm( sr1, N );
fr2  = NonlinearForcingTerm( sr2, N );

% Time scaling factor
tmp1 = 1.0; tmp2 = sqrt( 3 );

tscl1 = tau/tau_rd1 * tmp1;
tscl2 = tau/tau_rd2 * tmp2;

input_arr1 = fr1.calc_forcing_term( t_arr( 1:end-1 )/tscl1, rhy1.weight, t0i, eye( 2 ) );
input_arr2 = fr2.calc_forcing_term( t_arr( 1:end-1 )/tscl2, rhy2.weight, t0i, eye( 2 ) );

prim1 = input_arr1 + az*bz*g_rd1;
prim2 = input_arr2 + az*bz*g_rd2;

a = subplot( 3, 1, 1, 'parent', f );
hold on; axis equal

gain_arr = 0:0.1:1;
Ntmp = length( gain_arr );
off = 14.0;

scl1 =  0.4;
scl2 = 11.0;

for i = 1: Ntmp
    gain = gain_arr( Ntmp-i+1 );

    prim = gain * prim1 * scl1 + ( 1 - gain ) * prim2 * scl2;

    [ y_arr_comb, ~, ~] = ts.rollout( zeros( n, 1 ), zeros( n, 1 ), zeros( n, 1 ), prim, t0i, t_arr  );    

    plot( a, (i-1)*off+y_arr_comb( 1, : )-3, -(i-1)*0.020*off+y_arr_comb( 2, : )+3.5, 'linewidth', 3, 'color', c_orange )

end
set( a, 'xticklabel', {}, 'yticklabel', {}, 'xlim', [-15, 142], 'ylim', [-10, 10]  )
set( a, 'color', 'none' )


% Time scaling factor
tmp1 = 1.0; tmp2 = sqrt( 13 );

tscl1 = tau/tau_rd1 * tmp1;
tscl2 = tau/tau_rd2 * tmp2;

input_arr1 = fr1.calc_forcing_term( t_arr( 1:end-1 )/tscl1, rhy1.weight, t0i, eye( 2 ) );
input_arr2 = fr2.calc_forcing_term( t_arr( 1:end-1 )/tscl2, rhy2.weight, t0i, eye( 2 ) );

prim1 = input_arr1 + az*bz*g_rd1;
prim2 = input_arr2 + az*bz*g_rd2;

a = subplot( 3, 1, 2, 'parent', f );
hold on; axis equal

gain_arr = 0:0.1:1;
Ntmp = length( gain_arr );
off = 14.0;

scl1 =  0.4;
scl2 = 11.0;
tscl_arr = linspace( 1, 3, Ntmp );

for i = 1: Ntmp
    gain = gain_arr( Ntmp-i+1 );

    prim = gain * prim1 * scl1 + ( 1 - gain ) * prim2 * scl2;

    [ y_arr_comb, ~, ~] = ts.rollout( zeros( n, 1 ), zeros( n, 1 ), zeros( n, 1 ), prim, t0i, t_arr  );    

    plot( a, (i-1)*off+y_arr_comb( 1, : )-3, -(i-1)*0.020*off+y_arr_comb( 2, : )+3.5, 'linewidth', 3, 'color', c_orange )

end
set( a, 'xticklabel', {}, 'yticklabel', {}, 'xlim', [-15, 142], 'ylim', [-10, 10]  )
set( a, 'color', 'none' )


% Time scaling factor
tmp1 = 1.0; tmp2 = sqrt( 22 );

tscl1 = tau/tau_rd1 * tmp1;
tscl2 = tau/tau_rd2 * tmp2;

input_arr1 = fr1.calc_forcing_term( t_arr( 1:end-1 )/tscl1, rhy1.weight, t0i, eye( 2 ) );
input_arr2 = fr2.calc_forcing_term( t_arr( 1:end-1 )/tscl2, rhy2.weight, t0i, eye( 2 ) );

prim1 = input_arr1 + az*bz*g_rd1;
prim2 = input_arr2 + az*bz*g_rd2;

a = subplot( 3, 1, 3, 'parent', f );
hold on; axis equal

for i = 1: Ntmp
    gain = gain_arr( Ntmp-i+1 );

    prim = gain * prim1 * scl1 + ( 1 - gain ) * prim2 * scl2;

    [ y_arr_comb, ~, ~] = ts.rollout( zeros( n, 1 ), zeros( n, 1 ), zeros( n, 1 ), prim, t0i, t_arr  );    

    plot( a, (i-1)*off+y_arr_comb( 1, : )-3, -(i-1)*0.020*off+y_arr_comb( 2, : )+3.5, 'linewidth', 3, 'color', c_orange )

end
set( a, 'xticklabel', {}, 'yticklabel', {}, 'xlim', [-15, 142], 'ylim', [-10, 10]  )
set( a, 'color', 'none' )



fig_save( f, './images/fig2a_G' )



%% (1B) [Figure 2a] Sequencing Discrete Movements - Drawing Individual Images

close all; clc;

% The alphabets that we aim to draw
alphabets = { 'A', 'B', 'C' };
Na = length( alphabets );

% The trajectory data
traj_data = cell( 1, Na );

for i = 1 : Na
    a = alphabets{ i };

    % Call Data
    load( [ '../learned_parameters/discrete/', a, '_loose.mat' ] );
    traj_data{ i } = data;
end

traj_data{ 1 }.tau = traj_data{ 1 }.tau * 1.0;
traj_data{ 2 }.tau = traj_data{ 2 }.tau * 1.0;
traj_data{ 3 }.tau = traj_data{ 3 }.tau * 1.0;

Ntraj = Na;

% Define the time start 
t0_arr = traj_data{ 1 }.tau * (0:Na-1);

T     =   25.0;           % The   whole time of the simulation 
dt    =   1e-3;           % Time-step  for the simulation
t_arr = 0:dt:T;           % Time array for the simulation

% Position and velocity array of the trajectory
% This is 2D Data
 p_data_arr = zeros( 2, length( t_arr ), Ntraj );
dp_data_arr = zeros( 2, length( t_arr ), Ntraj );

% Scaling and Rotation of the trajectory
scl_arr = [ 1.3, 1.1, 1.1 ];

% The input force array for the movements
force_arr         = zeros( 2, length( t_arr )-1, Ntraj );
force_arr_trimmed = zeros( 2, length( t_arr )-1, Ntraj );

% The input force array for the movements
prim_for_plot = zeros( 2, length( t_arr )-1, Ntraj );

% The goal and initated time and final time for each movement
goal_arr = zeros( 2, Ntraj );

% Add position offset
xy_off = [ 0,0 ; 11.5, 11.5; 35, 12 ]'; 

% We need to first generate the DMP 
for i = 1 : Ntraj
    
    % Call the data from the trajectory
    data = traj_data{ i };
    
    % Getting the number of basis functions from the nonlinear forcing term
    [ ~, N ] = size( data.weight );

    % The Three elements of DMP
    cs        = CanonicalSystem( 'discrete', data.tau, data.alpha_s );
    fs        = NonlinearForcingTerm( cs, N );
    trans_sys = TransformationSystem( data.alpha_z, data.beta_z, cs );

    % Calculate the nonlinear forcing term 
    % This can be diminishing
    input_arr         = fs.calc_forcing_term( t_arr( 1:end-1 ), data.weight, t0_arr( i ), eye( 2 ) );
    input_arr_trimmed = fs.calc_forcing_term( t_arr( 1:end-1 ), data.weight, t0_arr( i ), eye( 2 ), 'trimmed' );

    % Scaling
    scl = scl_arr( i );
        
    % Calculate the new goal location and save
    new_goal = scl*data.goal;
    goal_arr( :, i ) = new_goal;

    % Calculate the new force and save
    new_force         = scl*input_arr; 
    new_force_trimmed = scl*input_arr_trimmed;

    % The nonlinear Forcing Terms
    force_arr( :, :, i )         = new_force;
    force_arr_trimmed( :, :, i ) = new_force_trimmed;

    prim_for_plot( :, :, i ) =  fs.calc_forcing_term( t_arr( 1:end-1 ), data.weight, 0, eye( 2 ) ) + data.alpha_z * data.beta_z * goal_arr( :, i );

    % Rollout, note that initial condition is set to be zeros. 
    [ y_arr, ~, dy_arr ] = trans_sys.rollout( zeros( 2, 1 ), zeros( 2, 1 ), new_goal, new_force, t0_arr( i ), t_arr  );    
    
    % These values are saved, just for sanity check 
     p_data_arr( :, :, i ) =  y_arr;
    dp_data_arr( :, :, i ) = dy_arr;

end

% The primitive forces 
prim         = zeros( 2, length( t_arr )-1, Ntraj );
prim_trimmed = zeros( 2, length( t_arr )-1, Ntraj );

% The primitives, including the goal offset.
for i = 1 : Ntraj
    prim( :, :, i )  = force_arr( :, :, i ) + data.alpha_z * data.beta_z * ( goal_arr( :, i ) + xy_off( :, i ) );
end


for i = 1 : Ntraj
f = figure( );
a1 = subplot( 1, 2, 1 );
hold on
plot( a1, p_data_arr( 1, :, i ), p_data_arr( 2, :, i ), 'linewidth', 10, 'color', c_blue )
scatter( a1, p_data_arr( 1,   1, i ), p_data_arr( 2,   1, i ), 500, 'filled', 'o', 'markerfacecolor', 'w', 'markeredgecolor', c_blue, 'linewidth', 4 )
scatter( a1, p_data_arr( 1, end, i ), p_data_arr( 2, end, i ), 500, 'filled', 'd', 'markerfacecolor', 'w', 'markeredgecolor', c_blue, 'linewidth', 4 )
axis equal
set( a1, 'xticklabel', {} , 'yticklabel', {} )
xlabel( '$X$ (m)', 'fontsize', 40 )
ylabel( '$Y$ (m)', 'fontsize', 40 )
tmpx = get( a1, 'xlim' ); tmpy = get( a1, 'ylim' );
tmpx( 1 ) = tmpx( 1 ) - 1.5; tmpx( 2 ) = tmpx( 2 ) + 1.5;
tmpy( 1 ) = tmpy( 1 ) - 1.5; tmpy( 2 ) = tmpy( 2 ) + 1.5;
set( a1, 'xlim', tmpx, 'ylim', tmpy )

a2 = subplot( 1, 2, 2 );
hold on
xlabel( '$t$ (s)', 'fontsize', 40 )
ylabel( '$\mathbf{f}_d(t)$ (-)', 'fontsize', 40 )
plot( t_arr( 1:end-1 ), prim_for_plot( 1, :, i ), 'linewidth', 10, 'color', c_blue, 'linestyle', '-' )
plot( t_arr( 1:end-1 ), prim_for_plot( 2, :, i ), 'linewidth', 10, 'color', c_blue, 'linestyle', ':' )
set( a2, 'xlim', [ 0, traj_data{ i }.tau ] )
set( a2, 'xticklabel', {} , 'yticklabel', {} )
legend( 'X', 'Y', 'fontsize', 40, 'location', 'northwest' )

fig_save( f, [ './images/fig2a_', alphabets{ i } ] )

end

%% (1C) [Figure 2b] Sequencing Discrete Movements - Merging Movements

f = figure( ); 
a1 = subplot( 3, 3, [1:3] );
hold on
% Defining the activation functions
act_arr = zeros( Ntraj, length( t_arr ) );
for i = 1 : Ntraj
    % If the first movement, 
    if i == 1
        t_end   = t0_arr( i )+traj_data{ i }.tau;
        act_arr( i, : ) = sigmoid_activation( t_arr, t_end-0.7, t_end+0.7, 1 );
    elseif i == 3
        t_start = t0_arr( i );
        act_arr( i, : ) = sigmoid_activation( t_arr, t_start-0.1, t_start+0.5, 0 );
    else
        t_start = t0_arr( i );
        t_end   = t0_arr( i )+traj_data{ i }.tau;
        act_arr( i, : ) =  sigmoid_activation( t_arr, t_start-0.5, t_start+0.5, 0 ) - sigmoid_activation( t_arr, t_end-0.3, t_end+0.5, 0 );
    end
end

act_arr = act_arr( :, 1:end-1 );

for i = 1 : Ntraj
    plot( a1, t_arr( 1:end-1 ), act_arr( i, : ), 'color', c_blue )
    area( a1, t_arr( 1:end-1 ), act_arr( i, : ), 'facecolor', c_blue, 'facealpha', 0.4, 'edgealpha', 0.0 )
end
set( a1, 'fontsize', 30, 'xlim', [ 0., 18.0 ] )
% xlabel( '$t$ (s)', 'fontsize', 40 )
% ylabel( '$\alpha_i(t)$(-)', 'fontsize', 40 )


a2 = subplot( 3, 3, [4:6] );
hold on
summed_prim = zeros( 2, length( t_arr ) - 1 );

for i = 1 : Ntraj
    summed_prim = summed_prim + prim( :, :, i ) .* act_arr( i, : );
end

plot( a2, t_arr( 1:end-1 ), summed_prim( 1, : ), 'linestyle', '-', 'color', c_blue )
plot( a2, t_arr( 1:end-1 ), summed_prim( 2, : ), 'linestyle', ':', 'color', c_blue )

set( a2, 'fontsize', 30, 'xlim', [ 0., 18.0 ], 'yticklabel', {}  )
xlabel( '$t$ (s)', 'fontsize', 40 )
% ylabel( '$\sum_{i=1}^{3}\alpha_i(t)\mathbf{f}_{d,i}(t)$(-)', 'fontsize', 40 )

tmp_arr = [ 5500, 11300, 25000 ];

a3 = subplot( 3, 3, 7 );
hold on

plot( a3, p_data_arr( 1, :, 1 ) + xy_off( 1, 1 ), p_data_arr( 2, :, 1 ) + xy_off( 2, 1 ), 'linewidth', 3, 'color', 0.3*ones(1,3), 'linestyle', ':' )
plot( a3, p_data_arr( 1, :, 2 ) + xy_off( 1, 2 ), p_data_arr( 2, :, 2 ) + xy_off( 2, 2 ), 'linewidth', 3, 'color', 0.3*ones(1,3), 'linestyle', ':' )
plot( a3, p_data_arr( 1, :, 3 ) + xy_off( 1, 3 ), p_data_arr( 2, :, 3 ) + xy_off( 2, 3 ), 'linewidth', 3, 'color', 0.3*ones(1,3), 'linestyle', ':' )

% The Three elements of DMP
cs        = CanonicalSystem( 'discrete', data.tau, data.alpha_s );
fs        = NonlinearForcingTerm( cs, N );
trans_sys = TransformationSystem( data.alpha_z, data.beta_z, cs );

% Rollout, note that initial condition is set to be zeros. 
[ y_arr, ~, ~] = trans_sys.rollout( zeros( 2, 1 ), zeros( 2, 1 ), zeros( 2, 1 ), summed_prim, 0.0, t_arr  );    
plot( a3, y_arr( 1, 1:tmp_arr( 1 ) ), y_arr( 2, 1:tmp_arr( 1 ) ), 'color', c_blue )
axis equal
ylabel( '$Y$ (m)', 'fontsize', 40 )
set( a3, 'xlim', [-3, 40], 'ylim', [-3, 18 ], 'xticklabel', {}, 'yticklabel', {} )

scatter( a3, p_data_arr( 1,   1, 1 ) + xy_off( 1, 1 ), p_data_arr( 2,   1, 1 ) + xy_off( 2, 1 ), 200,'filled', 'o', 'markerfacecolor', 'w', 'markeredgecolor', c_blue, 'linewidth', 5 )
% scatter( a3, p_data_arr( 1, end, 1 ) + xy_off( 1, 1 ), p_data_arr( 2, end, 1 ) + xy_off( 2, 1 ), 200,'filled', 'd', 'markerfacecolor', 'w', 'markeredgecolor', c_blue, 'linewidth', 5 )

% scatter( a3, p_data_arr( 1,   1, 2 ) + xy_off( 1, 2 ), p_data_arr( 2,   1, 2 ) + xy_off( 2, 2 ), 200, 'filled', 'o', 'markerfacecolor', 'w', 'markeredgecolor', c_blue, 'linewidth', 5 )
% scatter( a3, p_data_arr( 1, end, 2 ) + xy_off( 1, 2 ), p_data_arr( 2, end, 2 ) + xy_off( 2, 2 ), 200,'filled', 'd', 'markerfacecolor', 'w', 'markeredgecolor', c_blue, 'linewidth', 5 )
% 
% scatter( a3, p_data_arr( 1,   1, 3 ) + xy_off( 1, 3 ), p_data_arr( 2,   1, 3 ) + xy_off( 2, 3 ), 200,'filled', 'o', 'markerfacecolor', 'w', 'markeredgecolor', c_blue, 'linewidth', 5 )
% scatter( a3, p_data_arr( 1, end, 3 ) + xy_off( 1, 3 ), p_data_arr( 2, end, 3 ) + xy_off( 2, 3 ), 200,'filled', 'd', 'markerfacecolor', 'w', 'markeredgecolor', c_blue, 'linewidth', 5 )

a3 = subplot( 3, 3, 8 );
hold on

plot( a3, p_data_arr( 1, :, 1 ) + xy_off( 1, 1 ), p_data_arr( 2, :, 1 ) + xy_off( 2, 1 ), 'linewidth', 3, 'color', 0.3*ones(1,3), 'linestyle', ':' )
plot( a3, p_data_arr( 1, :, 2 ) + xy_off( 1, 2 ), p_data_arr( 2, :, 2 ) + xy_off( 2, 2 ), 'linewidth', 3, 'color', 0.3*ones(1,3), 'linestyle', ':' )
plot( a3, p_data_arr( 1, :, 3 ) + xy_off( 1, 3 ), p_data_arr( 2, :, 3 ) + xy_off( 2, 3 ), 'linewidth', 3, 'color', 0.3*ones(1,3), 'linestyle', ':' )


% The Three elements of DMP
cs        = CanonicalSystem( 'discrete', data.tau, data.alpha_s );
fs        = NonlinearForcingTerm( cs, N );
trans_sys = TransformationSystem( data.alpha_z, data.beta_z, cs );

% Rollout, note that initial condition is set to be zeros. 
[ y_arr, ~, ~ ] = trans_sys.rollout( zeros( 2, 1 ), zeros( 2, 1 ), zeros( 2, 1 ), summed_prim, 0.0, t_arr  );    

plot( a3, y_arr( 1, 1:tmp_arr( 2 ) ), y_arr( 2, 1:tmp_arr( 2 )  ), 'color', c_blue )
axis equal
xlabel( '$X$ (m)', 'fontsize', 40 )
set( a3, 'xlim', [-3, 40], 'ylim', [-3, 18 ], 'xticklabel', {}, 'yticklabel', {} )


scatter( a3, p_data_arr( 1,   1, 1 ) + xy_off( 1, 1 ), p_data_arr( 2,   1, 1 ) + xy_off( 2, 1 ), 200,'filled', 'o', 'markerfacecolor', 'w', 'markeredgecolor', c_blue, 'linewidth', 5 )
% scatter( a3, p_data_arr( 1, end, 1 ) + xy_off( 1, 1 ), p_data_arr( 2, end, 1 ) + xy_off( 2, 1 ), 200,'filled', 'd', 'markerfacecolor', 'w', 'markeredgecolor', c_blue, 'linewidth', 5 )

% scatter( a3, p_data_arr( 1,   1, 2 ) + xy_off( 1, 2 ), p_data_arr( 2,   1, 2 ) + xy_off( 2, 2 ), 200, 'filled', 'o', 'markerfacecolor', 'w', 'markeredgecolor', c_blue, 'linewidth', 5 )
% scatter( a3, p_data_arr( 1, end, 2 ) + xy_off( 1, 2 ), p_data_arr( 2, end, 2 ) + xy_off( 2, 2 ), 200,'filled', 'd', 'markerfacecolor', 'w', 'markeredgecolor', c_blue, 'linewidth', 5 )

% scatter( a3, p_data_arr( 1,   1, 3 ) + xy_off( 1, 3 ), p_data_arr( 2,   1, 3 ) + xy_off( 2, 3 ), 200,'filled', 'o', 'markerfacecolor', 'w', 'markeredgecolor', c_blue, 'linewidth', 5 )
% scatter( a3, p_data_arr( 1, end, 3 ) + xy_off( 1, 3 ), p_data_arr( 2, end, 3 ) + xy_off( 2, 3 ), 200,'filled', 'd', 'markerfacecolor', 'w', 'markeredgecolor', c_blue, 'linewidth', 5 )

a3 = subplot( 3, 3, 9 );
hold on

plot( a3, p_data_arr( 1, :, 1 ) + xy_off( 1, 1 ), p_data_arr( 2, :, 1 ) + xy_off( 2, 1 ), 'linewidth', 3, 'color', 0.3*ones(1,3), 'linestyle', ':' )
plot( a3, p_data_arr( 1, :, 2 ) + xy_off( 1, 2 ), p_data_arr( 2, :, 2 ) + xy_off( 2, 2 ), 'linewidth', 3, 'color', 0.3*ones(1,3), 'linestyle', ':' )
plot( a3, p_data_arr( 1, :, 3 ) + xy_off( 1, 3 ), p_data_arr( 2, :, 3 ) + xy_off( 2, 3 ), 'linewidth', 3, 'color', 0.3*ones(1,3), 'linestyle', ':' )

% The Three elements of DMP
cs        = CanonicalSystem( 'discrete', data.tau, data.alpha_s );
fs        = NonlinearForcingTerm( cs, N );
trans_sys = TransformationSystem( data.alpha_z, data.beta_z, cs );

% Rollout, note that initial condition is set to be zeros. 
[ y_arr, ~, ~] = trans_sys.rollout( zeros( 2, 1 ), zeros( 2, 1 ), zeros( 2, 1 ), summed_prim, 0.0, t_arr  );    

plot( a3, y_arr( 1, 1:tmp_arr( 3 )  ), y_arr( 2, 1:tmp_arr( 3 ) ), 'color', c_blue )
axis equal
set( a3, 'xlim', [-3, 40], 'ylim', [-3, 18 ], 'xticklabel', {}, 'yticklabel', {} )


scatter( a3, p_data_arr( 1,   1, 1 ) + xy_off( 1, 1 ), p_data_arr( 2,   1, 1 ) + xy_off( 2, 1 ), 200,'filled', 'o', 'markerfacecolor', 'w', 'markeredgecolor', c_blue, 'linewidth', 5 )
% scatter( a3, p_data_arr( 1, end, 1 ) + xy_off( 1, 1 ), p_data_arr( 2, end, 1 ) + xy_off( 2, 1 ), 200,'filled', 'd', 'markerfacecolor', 'w', 'markeredgecolor', c_blue, 'linewidth', 5 )
% 
% scatter( a3, p_data_arr( 1,   1, 2 ) + xy_off( 1, 2 ), p_data_arr( 2,   1, 2 ) + xy_off( 2, 2 ), 200, 'filled', 'o', 'markerfacecolor', 'w', 'markeredgecolor', c_blue, 'linewidth', 5 )
% scatter( a3, p_data_arr( 1, end, 2 ) + xy_off( 1, 2 ), p_data_arr( 2, end, 2 ) + xy_off( 2, 2 ), 200,'filled', 'd', 'markerfacecolor', 'w', 'markeredgecolor', c_blue, 'linewidth', 5 )

% scatter( a3, p_data_arr( 1,   1, 3 ) + xy_off( 1, 3 ), p_data_arr( 2,   1, 3 ) + xy_off( 2, 3 ), 200,'filled', 'o', 'markerfacecolor', 'w', 'markeredgecolor', c_blue, 'linewidth', 5 )
scatter( a3, p_data_arr( 1, end, 3 ) + xy_off( 1, 3 ), p_data_arr( 2, end, 3 ) + xy_off( 2, 3 ), 200,'filled', 'd', 'markerfacecolor', 'w', 'markeredgecolor', c_blue, 'linewidth', 5 )

fig_save( f, './images/fig2b')


%% (1D) [Figure 3a] Merging the Primitives

% Merging Two Discrete Movements
close all; clc;

f = figure( ); 

% Load two mat functions, discrete movements
tmp1 = load( '../learned_parameters/discrete/A_loose.mat' ); dis1 = tmp1.data;
tmp2 = load( '../learned_parameters/discrete/B_loose.mat' ); dis2 = tmp2.data;

% Get the degrees of freedom and the number of weights 
n = size( dis1.weight, 1 );
N = size( dis1.weight, 2 );

% The alpha, beta values must match
% The Three elements of DMP
cs        = CanonicalSystem( 'discrete', dis1.tau, dis1.alpha_s );
fs        = NonlinearForcingTerm( cs, N );
trans_sys = TransformationSystem( dis1.alpha_z, dis1.beta_z, cs );

% Rollout for each discrete movement
t0i   =    0.0;           % The initial time of the movement rollout
dt    =   1e-4;           % Time-step  for the simulation
t_arr = 0:dt:(2*pi*3);    % Time array for the simulation

F_arr1 = fs.calc_forcing_term( t_arr( 1:end-1 ), dis1.weight, t0i, eye( 2 ) );
F_arr2 = fs.calc_forcing_term( t_arr( 1:end-1 ), dis2.weight, t0i, eye( 2 ) );
 
a = subplot( 3, 1, 1, 'parent', f );
hold on; axis equal

gain_arr = 0:0.1:1;
Ntmp = length( gain_arr );
off = 12.0;

scl1 = 1.1;
scl2 = 0.9;

for i = 1: Ntmp
    gain = gain_arr( Ntmp-i+1 );

    F_arr =  scl1 * F_arr1 * gain    +  scl2* F_arr2 * ( 1 - gain );
    goal  =  scl1 * dis1.goal * gain + scl2* dis2.goal * ( 1 - gain );

    [ y_arr_comb, ~, ~] = trans_sys.rollout( zeros( n, 1 ), zeros( n, 1 ), zeros( n, 1 ), F_arr + dis1.alpha_z*dis1.beta_z*goal, t0i, t_arr  );    
    
    plot( a, (i-1)*off+y_arr_comb( 1, : )-3, (i-1)*0.08*off+y_arr_comb( 2, : )+3.5, 'linewidth', 6, 'color', c_blue )
end
set( a, 'xticklabel', {}, 'yticklabel', {}, 'xlim', [-12, 132], 'ylim', [0, 20]  )
set( a, 'color', 'none' )

% Load two mat functions, discrete movements
tmp1 = load( '../learned_parameters/discrete/A_loose.mat' ); dis1 = tmp1.data;
tmp2 = load( '../learned_parameters/rhythmic/heart.mat'   ); rhy1 = tmp2.data;

% Get the degrees of freedom and the number of weights 
n = size( dis1.weight, 1 );
N = size( dis2.weight, 2 );

% Discrete DMP
dis1.tau = 1.0;
cs_d        = CanonicalSystem( 'discrete', dis1.tau, dis1.alpha_s );
fs_d        = NonlinearForcingTerm( cs_d, N );
trans_sys_d = TransformationSystem( dis1.alpha_z, dis1.beta_z, cs_d );

% Rhythmic DMP
cs_r        = CanonicalSystem( 'rhythmic', rhy1.tau, dis1.alpha_s );
fs_r        = NonlinearForcingTerm( cs_r, N );
trans_sys_r = TransformationSystem( rhy1.alpha_z, rhy1.beta_z, cs_r );

% Rollout for each discrete movement
t0i   =    0.0;           % The initial time of the movement rollout
T     =  dis1.tau;        % The   whole time of the simulation 
dt    =   1e-3;           % Time-step  for the simulation
t_arr = 0:dt:(1*T);       % Time array for the simulation

F_arr1 = fs_d.calc_forcing_term( t_arr( 1:end-1 ), dis1.weight, t0i, eye( 2 ) );
F_arr2 = fs_r.calc_forcing_term( t_arr( 1:end-1 ), rhy1.weight, t0i, eye( 2 ) );
 
a = subplot( 3, 1, 2, 'parent', f );
hold on; axis equal

gain_arr = 0:0.1:1;
Ntmp = length( gain_arr );
off = 12.5;

% Initialize a matrix to hold the interpolated colors
color_arr = zeros( Ntmp, 3);

% Interpolate each color channel
for i = 1:3
    color_arr(:, i) = linspace( c_orange( i ), c_blue( i ), Ntmp );
end

scl1 = 1.1;
scl2 = 0.5;

for i = 1: Ntmp
    gain = gain_arr( Ntmp-i+1 );

    F_arr =    scl1*F_arr1 * gain +   scl2*F_arr2 * ( 1 - gain );
    goal  = scl1*dis1.goal * gain +   scl2*( rhy1.goal - rhy1.p_init )* ( 1 - gain );
   
    [ y_arr_comb, ~, ~] = trans_sys_d.rollout( zeros( 2, 1 ), zeros( 2, 1), zeros( 2, 1), F_arr + dis1.alpha_z*dis1.beta_z*goal, t0i, t_arr  );  
    
    plot( a, (i-1)*off+y_arr_comb( 1, : )-3, (i-1)*0.08*off+y_arr_comb( 2, : )+3.5, 'linewidth', 6, 'color', color_arr( Ntmp-i+1, : ) )

end
set( a, 'xticklabel', {}, 'yticklabel', {}, 'xlim', [-12, 132], 'ylim', [0, 20]  )
set( a, 'color', 'none' )

% Load two mat functions, discrete movements
tmp1 = load( '../learned_parameters/rhythmic/circle.mat' ); rhy1 = tmp1.data;
tmp2 = load( '../learned_parameters/rhythmic/heart.mat'  ); rhy2 = tmp2.data;

% Get the degrees of freedom and the number of weights 
n = size( rhy1.weight, 1 );
N = size( rhy1.weight, 2 );

% Discrete DMP
cs_r1       = CanonicalSystem( 'rhythmic', rhy1.tau, 1. );
fs_r1       = NonlinearForcingTerm( cs_r1, N );
trans_sys_r1 = TransformationSystem( rhy1.alpha_z, rhy1.beta_z, cs_r1 );

% Rhythmic DMP
cs_r2        = CanonicalSystem( 'rhythmic', rhy2.tau, 1. );
fs_r2        = NonlinearForcingTerm( cs_r2, N );
trans_sys_r2 = TransformationSystem( rhy2.alpha_z, rhy2.beta_z, cs_r2 );

% Rollout for each discrete movement
t0i   =    0.0;           % The initial time of the movement rollout
T     =   2*pi*6;           % The   whole time of the simulation 
dt    =   1e-3;           % Time-step  for the simulation
t_arr = 0:dt:T;           % Time array for the simulation

scl1 =  3.0;
scl2 =  0.6;

F_arr1 = fs_r1.calc_forcing_term( t_arr( 1:end-1 ), rhy1.weight, t0i, scl1*eye( 2 ) );
F_arr2 = fs_r2.calc_forcing_term( t_arr( 1:end-1 ), rhy2.weight, t0i, scl2*eye( 2 ) );
 
a = subplot( 3, 1, 3, 'parent', f );
hold on; axis equal

gain_arr = 0:0.1:1;
Ntmp = length( gain_arr );
off = 12.0;

% Initialize a matrix to hold the interpolated colors
color_arr = zeros( Ntmp, 3);

for i = 1: Ntmp
    gain = gain_arr( i );

    F_arr =    F_arr1 * gain * scl1 +   F_arr2 * ( 1 - gain )* scl2;
    goal  = rhy1.goal * gain * scl1+   rhy2.goal * ( 1 - gain )* scl2;
     init =  rhy1.p_init * gain * scl1+ rhy2.p_init * ( 1 - gain )* scl2;
    zinit =  rhy1.dp_init*rhy1.tau * gain* scl1 +  rhy2.dp_init*rhy2.tau * ( 1-gain ) * scl2;
   
    [ y_arr_comb, ~, ~] = trans_sys_r1.rollout( init, zinit, goal, F_arr, t0i, t_arr  );    
    plot( a, (i-1)*off+y_arr_comb( 1, 500:end ), -(i-1)*0.00*off+y_arr_comb( 2, 500:end ), 'linewidth', 6, 'color', c_orange )

end
set( a, 'xticklabel', {}, 'yticklabel', {}, 'xlim', [-12, 132], 'ylim', [-10, 10]  )
set( a, 'color', 'none' )
fig_save( f, './images/fig3c')

%% (1D) [Figure 3a] Merging the Primitives, Chaotic

% Merging Two Rhythmic Movements

% Load two mat functions, discrete movements
tmp1 = load( '../learned_parameters/rhythmic/heart.mat' );  rhy1 = tmp1.data;
tmp2 = load( '../learned_parameters/rhythmic/circle.mat' ); rhy2 = tmp2.data;

w_ratio = [ sqrt( 3 ), sqrt( 17 ), sqrt( 22 ) ];
tmp = { 'a', 'b', 'c' }; 
my_title = { '$\tau_1/\tau_2=\sqrt{3}$', '$\tau_1/\tau_2=\sqrt{17}$', '$\tau_1/\tau_2=\sqrt{22}$' };
f = figure( );
for j = 1:length( w_ratio )

a = subplot( length( w_ratio ), 1, j )
hold on

% Get the degrees of freedom and the number of weights 
n = size( rhy1.weight, 1 );
N = size( rhy1.weight, 2 );

% Create DMP for rhythmic
cs_r1        = CanonicalSystem( 'rhythmic', rhy1.tau*w_ratio( j ), 1.0 );
fs_r1        = NonlinearForcingTerm( cs_r1, N );
trans_sys_r1 = TransformationSystem( rhy1.alpha_z, rhy1.beta_z, cs_r1 );


% Create DMP for rhythmic
cs_r2        = CanonicalSystem( 'rhythmic', rhy2.tau, 1.0 );
fs_r2        = NonlinearForcingTerm( cs_r2, N );
trans_sys_r2 = TransformationSystem( rhy2.alpha_z, rhy2.beta_z, cs_r1 );

% Rollout for each discrete movement
T     =   2*pi*4;           % The   whole time of the simulation 
dt    =   1e-3;           % Time-step  for the simulation
t_arr = 0:dt:T;           % Time array for the simulation

scl1 = 1;
scl2 = 26;

F_arr1 = fs_r1.calc_forcing_term( t_arr( 1:end-1 ), rhy1.weight, 0, scl1 * eye( 2 ) );
F_arr2 = fs_r2.calc_forcing_term( t_arr( 1:end-1 ), rhy2.weight, 0, scl2 * eye( 2 ) );

gain_arr = 0:0.1:1;
Ntmp = length( gain_arr );
off = 33.0;

for i = 1: Ntmp
    gain = gain_arr( i );

    F_arr =    F_arr1 * gain +    F_arr2 * ( 1 - gain );
     init = scl1 *  rhy1.p_init * gain + scl2 * rhy2.p_init * (1 - gain );
    dinit = scl1 * rhy1.dp_init * gain + scl2 * rhy2.dp_init * (1 - gain );
    
    
    [ y_arr_comb, ~, ~] = trans_sys_r1.rollout( init, dinit, zeros( 2, 1), F_arr, 0, t_arr  );    
    
    plot( a, (i-1)*off+y_arr_comb( 1, : ), y_arr_comb( 2, : ), 'linewidth', 2, 'color', c_orange )
end
axis equal
set( gca,'xticklabel', {}, 'yticklabel', {}, 'xlim', [-25, 360], 'ylim', [-40, 40] )
title( my_title{ j }, 'fontsize', 40 )

end

fig_save( f, './images/fig4' )
