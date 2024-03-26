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
% Part1: Discrete DMP
% ==================================== %
tmp  = load( '../learned_parameters/discrete/A_loose.mat' );
data = tmp.data;

% The number of basis functions N
[ ~, N ] = size( data.weight );

% The three elements of discrete DMP
cs1        = CanonicalSystem( 'discrete', data.tau, data.alpha_s );
fs1        = NonlinearForcingTerm( cs1, N );
trans_sys1 = TransformationSystem( data.alpha_z, data.beta_z, cs1 );

% The Parameters for Forward Simulation
t0i   = 0.0;
T     = 16;
dt    = 1e-2;
t_arr = 0:dt:T;

% Calculate the nonlinear forcing term for discrete movement and rollout
input_arr = fs1.calc_forcing_term( t_arr( 1:end-1 ), data.weight, t0i, eye( 2 ) );
[ y_arr, ~, ~ ] = trans_sys1.rollout( zeros( 2, 1 ), zeros( 2, 1 ), zeros( 2, 1 ), data.alpha_z*data.beta_z*data.goal+input_arr, t0i, t_arr  );  

% Drawing the plots
f = figure( ); 
a1 = subplot( 3, 2, 1 );
plot( a1, t_arr, cs1.calc( t_arr ), 'linewidth', 4, 'color',  c_blue );
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

% For rhythmic movement, load the data
tmp  = load( '../learned_parameters/rhythmic/heart.mat' );
data = tmp.data;

[ ~, N ] = size( data.weight );

% The three elements of rhythmic DMP
% Note that this definition of the Canonical System is the simple linear function, 
% as we don't have a network of rhythmic Canonical System
cs2        = CanonicalSystem( 'rhythmic', data.tau, 1.0 );
fs2        = NonlinearForcingTerm( cs2, N );
trans_sys2 = TransformationSystem( data.alpha_z, data.beta_z, cs2 );

% The parameters for forward simulation
t0i   = 0.0;
T     = 3;
dt    = 1e-4;
t_arr = 0:dt:T;
Nt    = length( t_arr );

% Calculate the nonlinear forcing term for rhythmic movement and rollout
input_arr2 = fs2.calc_forcing_term( t_arr( 1:end-1 ), data.weight , t0i, eye( 2 ) );
off = [ 0.0; 0.8 ];
[ y_arr2, ~, ~ ] = trans_sys2.rollout( off , data.dp_init*data.tau, zeros( 2, 1 ), data.alpha_z*data.beta_z*(data.goal+off)+input_arr2, t0i, t_arr  );

a2 = subplot( 3, 2, 2 );
plot( t_arr, cs2.calc( t_arr ), 'linewidth', 4, 'color', c_orange);
set( a2,  'xlim', [0,3], 'xtick', 0:1:3, 'fontsize', 30, 'ylim', [ 0, 2*pi], 'ytick', 2*pi*[0., 0.5, 1.0], 'yticklabel', {'0', '$\pi$', '$2\pi$'} )
ylabel( '$s_r(t)$', 'fontsize', 40 )

a4 = subplot( 3, 2, 4 );
hold on
plot( a4, t_arr( 1:end-1 ), input_arr2( 1, : ), 'linewidth', 4, 'linestyle', '-' , 'color', c_orange );
plot( a4, t_arr( 1:end-1 ), input_arr2( 2, : ), 'linewidth', 4, 'linestyle', '-' , 'color', c_orange );
set( a4, 'yticklabel', {}, 'xlim', [0,3], 'xtick', 0:1:3, 'fontsize', 30 )
xlabel( '$t$', 'fontsize', 40 )
ylabel( '$\mathbf{f}_r(t)$', 'fontsize', 40 )

a6 = subplot( 3, 2, 6 );
hold on
plot( a6, y_arr2( 1, 1:end-1 ), y_arr2( 2, 1:end-1 ), 'linewidth', 4, 'color', c_orange );
scatter( a6, 0, 0, 100, 'filled',  'o', 'markerfacecolor', 'w', 'markeredgecolor', c_orange, 'linewidth', 4 )

set( a6, 'xticklabel', {}, 'yticklabel', {}, 'xlim', [-20,20], 'ylim',[-25,15])

xlabel( '$X$', 'fontsize', 40 )
ylabel( '$Y$', 'fontsize', 40 )
axis equal

fig_save( f, './images/fig1' )

%% (1Aa) [Figure 1a] Spatial Scaling

clear data*, close all; clc;
% Load the A alphabet 
tmp  = load( '../learned_parameters/discrete/A_loose.mat' );
data_d = tmp.data;

% The number of basis functions N
[ ~, N ] = size( data_d.weight );

% The three elements of discrete DMP
cs_d        = CanonicalSystem( 'discrete', data_d.tau, data_d.alpha_s );
fs_d        = NonlinearForcingTerm( cs_d, N );
trans_sys_d = TransformationSystem( data_d.alpha_z, data_d.beta_z, cs_d );

f = figure( ); 
a1 = subplot( 1, 2, 1 );
hold on; axis equal
scl_arr = [ 0.3, 0.5, 0.7, 0.9, 1.0, 1.2, 1.5, 2.0 ];
offset_arr = 10 * [ 0.2, 0.46, 0.9, 1.45, 2.2, 3.0, 4.0, 5.2];

% The Parameters for Forward Simulation
t0i   = 0.0;
T     = data_d.tau*2;
dt    = 1e-4;
t_arr = 0:dt:T;

for i = 1 : length( scl_arr )
    scl = scl_arr( i );
    off = offset_arr( i );

    % Calculate the nonlinear forcing term for discrete movement and rollout
    input_arr = fs_d.calc_forcing_term( t_arr( 1:end-1 ), data_d.weight, t0i, eye( 2 ) );
   
    % initial position
    g  = scl*( data_d.goal );
    % Traditional
    % [ y_arr, ~, ~ ] = trans_sys1.rollout( y0, data.z0, g,scl*input_arr, t0i, t_arr  );  
    [ y_arr, ~, ~ ] = trans_sys_d.rollout( zeros( 2, 1 ), zeros( 2, 1 ), zeros( 2, 1 ), scl*input_arr +data_d.alpha_z*data_d.beta_z*g, t0i, t_arr  );  
    % New one
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

a2 = subplot( 1, 2, 2 );
hold on; axis equal

% For rhythmic movement, load the data
tmp  = load( '../learned_parameters/rhythmic/heart.mat' );
data_r = tmp.data;

[ ~, N ] = size( data_r.weight );

% The three elements of rhythmic DMP
% Note that this definition of the Canonical System is the simple linear function, 
% as we don't have a network of rhythmic Canonical System
cs_r        = CanonicalSystem( 'rhythmic', data_r.tau, 1.0 );
fs_r        = NonlinearForcingTerm( cs_r, N );
trans_sys_r = TransformationSystem( data_r.alpha_z, data_r.beta_z, cs_r );

% The Parameters for Forward Simulation
t0i   = 0.0;
T     = data_r.tau*2*pi;
dt    = 1e-4;
t_arr = 0:dt:T;

off = 0;

for i = 1 : length( scl_arr )
    scl = scl_arr( i );

    % Calculate the nonlinear forcing term for discrete movement and rollout
    input_arr = fs_r.calc_forcing_term( t_arr( 1:end-1 ), data_r.weight, t0i, eye( 2 ) );
   
    % initial position
    g  = scl*data_r.goal ;
    % Traditional
    % [ y_arr, ~, ~ ] = trans_sys1.rollout( y0, data.z0, g,scl*input_arr, t0i, t_arr  );  
    [ y_arr, ~, ~ ] = trans_sys_r.rollout( y0, zeros( 2, 1 ), zeros( 2, 1 ), scl*input_arr+data_r.alpha_z*data_r.beta_z*g, t0i, t_arr  );  

    if scl ~= 1
        plot( a2, off + y_arr( 1, : ), y_arr( 2,  : ), 'color', c_orange, 'linewidth', 3 )
    else
        plot( a2, off + y_arr( 1, : ), y_arr( 2,  : ), 'color', c_orange, 'linewidth', 8 )
    end
end

scatter( a2, 0, 0, 100, 'filled',  'o', 'markerfacecolor', 'w', 'markeredgecolor', c_orange, 'linewidth', 3 )
xlabel( 'X', 'fontsize', 40 )
ylabel( 'Y', 'fontsize', 40 )

set( a2, 'xlim', [-36, 36 ], 'ylim', [-51, 21], 'xticklabel', {}, 'yticklabel', {} )
fig_save( f, './images/fig1a' )

%% (1Aa) [Figure 1b] Rotational Scaling

clear data*, close all; clc;
% Load the A alphabet 
tmp  = load( '../learned_parameters/discrete/A_loose.mat' );
data_d = tmp.data;

% The number of basis functions N
[ ~, N ] = size( data_d.weight );

% The three elements of discrete DMP
cs_d        = CanonicalSystem( 'discrete', data_d.tau, data_d.alpha_s );
fs_d        = NonlinearForcingTerm( cs_d, N );
trans_sys_d = TransformationSystem( data_d.alpha_z, data_d.beta_z, cs_d );

f = figure( ); 
a1 = subplot( 1, 2, 1 );
hold on; axis equal
ang_arr = linspace( 0, 2*pi, 6 );
ang_arr = ang_arr( 1:end-1 );

% The Parameters for Forward Simulation
t0i   = 0.0;
T     = data_d.tau;
dt    = 1e-4;
t_arr = 0:dt:T;
scl = 2.5;
for i = 1 : length( ang_arr )

    ang = ang_arr( i );

    % Calculate the nonlinear forcing term for discrete movement and rollout
    input_arr = fs_d.calc_forcing_term( t_arr( 1:end-1 ), data_d.weight, t0i, eye( 2 ) );
 
    % Rotation matrices 
    R = [ cos( ang ), -sin( ang );
          sin( ang ),  cos( ang )];

    % initial position
    y0 = R*[10;0];
    g  = scl*R*( data_d.goal );

    [ y_arr, ~, ~ ] = trans_sys_d.rollout( zeros( 2, 1 ), zeros( 2, 1 ), zeros( 2, 1 ), R*scl*input_arr + data_d.alpha_z*data_d.beta_z*g, t0i, t_arr  );  
    % New one

    plot( a1, y0( 1 ) + y_arr( 1, : ), y0( 2 ) + y_arr( 2,  : ), 'color', c_blue )

    scatter( a1, y0( 1 ) +y_arr( 1,   1 ),y0( 2 ) + y_arr( 2,   1 ), 200, 'filled',  'o', 'markerfacecolor', 'w', 'markeredgecolor', c_blue, 'linewidth', 3 )
    scatter( a1, y0( 1 ) +y_arr( 1, end ), y0( 2 ) +y_arr( 2, end ), 200, 'filled',  'd', 'markerfacecolor', 'w', 'markeredgecolor', c_blue, 'linewidth', 3 )

end
set( a1, 'xlim', get( a1, 'xlim' ) + [-4, 2], 'ylim', get( a1, 'ylim' ) + [-3, 3], 'xticklabel',{},'yticklabel',{} )
xlabel( a1, 'X', 'fontsize', 40)
ylabel( a1, 'Y', 'fontsize', 40)

a2 = subplot( 1, 2, 2 );
hold on; axis equal

% For rhythmic movement, load the data
tmp  = load( '../learned_parameters/rhythmic/heart.mat' );
data_r = tmp.data;

[ ~, N ] = size( data_r.weight );

% The three elements of rhythmic DMP
% Note that this definition of the Canonical System is the simple linear function, 
% as we don't have a network of rhythmic Canonical System
cs_r        = CanonicalSystem( 'rhythmic', data_r.tau, 1.0 );
fs_r        = NonlinearForcingTerm( cs_r, N );
trans_sys_r = TransformationSystem( data_r.alpha_z, data_r.beta_z, cs_r );

% The Parameters for Forward Simulation
t0i   = 0.0;
T     = data_r.tau*2*pi;
dt    = 1e-4;
t_arr = 0:dt:T;
scl = 0.6;
for i = 1 : length( ang_arr )
    ang = ang_arr( i );

    % Rotation matrices 
    R = [ cos( ang ), -sin( ang );
          sin( ang ),  cos( ang )];

    % initial position
    y0 = R*[15;0];
    g  = scl*R*data_r.goal;

    % Calculate the nonlinear forcing term for discrete movement and rollout
    input_arr = fs_r.calc_forcing_term( t_arr( 1:end-1 ), data_r.weight, t0i, eye( 2 ) );
 
    [ y_arr, ~, ~ ] = trans_sys_r.rollout( zeros( 2, 1 ), zeros( 2, 1 ), zeros( 2, 1 ), R*scl*input_arr + data_r.alpha_z*data_r.beta_z*g, t0i, t_arr  );  
    % New one

    plot( a2, y0( 1 ) + y_arr( 1, : ), y0( 2 ) + y_arr( 2,  : ), 'color', c_orange )
    scatter( a2, y0( 1 ),y0( 2 ), 200, 'filled',  'o', 'markerfacecolor', 'w', 'markeredgecolor', c_orange, 'linewidth', 3 )
    
end



set( a2, 'xlim', get( a2, 'xlim' ) + [-4, 2], 'ylim', get( a2, 'ylim' ) + [-3, 3], 'xticklabel',{},'yticklabel',{} )
xlabel( a2, 'X', 'fontsize', 40)
ylabel( a2, 'Y', 'fontsize', 40)

fig_save( f, './images/fig1b' )

%% (1Ac) [Figure 1c] Temporal Scaling 

clear data*, close all; clc;
f = figure( );
% Load the A alphabet 
tmp  = load( '../learned_parameters/discrete/A_loose.mat' );
data_d = tmp.data;

% The number of basis functions N
[ ~, N ] = size( data_d.weight );

% The three elements of discrete DMP
cs_d = CanonicalSystem( 'discrete', data_d.tau, data_d.alpha_s );
fs_d = NonlinearForcingTerm( cs_d, N );

% First, generate the trajectories
tscl_arr = [ 0.5, 1.0, 2.0 ];

% The Parameters for Forward Simulation
t0i   = 0.0;
T     = data_d.tau;
dt    = 1e-3;
t_arr = 0:dt:(T*2);

y_data_arr = zeros( 2, length( t_arr ), length( tscl_arr ) );

for i = 1 : length( tscl_arr )
    tscl = tscl_arr( i );

    % Rollout for data
    cs_tmp        = CanonicalSystem( 'discrete', tscl*data_d.tau, data_d.alpha_s );
    fs_tmp        = NonlinearForcingTerm( cs_tmp, N );
    trans_sys_tmp = TransformationSystem( data_d.alpha_z, data_d.beta_z, cs_tmp );
        
    % Calculate the nonlinear forcing term for discrete movement and rollout
    input_arr = fs_tmp.calc_forcing_term( t_arr( 1:end-1 ), data_d.weight, t0i, eye( 2 ) );
  
    [ y_arr, ~, ~ ] = trans_sys_tmp.rollout( zeros( 2, 1 ), zeros( 2, 1 ), zeros( 2, 1), input_arr +data_d.alpha_z*data_d.beta_z*data_d.goal, t0i, t_arr  );  

    y_data_arr( :, :, i ) = y_arr;
    
end

ttmp_arr = data_d.tau *(0.25:0.25:2.00);
offset  = 10*( 1:length( ttmp_arr ) );

% Going through the iterations
for i = 1: length( tscl_arr )
    a = subplot( length( tscl_arr ), 2, 2*i-1, 'parent', f );
    hold on; axis equal
    y_tmp = y_data_arr( :, :, i );
    for j = 1 : length( ttmp_arr )
        tt = ttmp_arr( j );

        [ ~, idx] = min( abs( t_arr - tt ) );
        off = offset( j );

        plot( a, off + y_tmp( 1, 1:idx ), y_tmp( 2, 1:idx ), 'linewidth', 5, 'color', c_blue )
        scatter( a, off + y_tmp( 1, idx ), y_tmp( 2, idx ), 100, 'markerfacecolor', 'w', 'markeredgecolor', c_blue, 'linewidth', 3 )        
    end
    set( a, 'xlim', [3, max(offset)+10], 'xticklabel', {}, 'yticklabel', {}  )
    
    title( a, ['$\tau=', num2str( tscl_arr( i ), '%.2f' ), '\tau^{(d)}$'], 'fontsize', 40 )
end

% Load the A alphabet 
tmp  = load( '../learned_parameters/rhythmic/heart.mat' );
data_r = tmp.data;

% The number of basis functions N
[ ~, N ] = size( data_r.weight );

% The three elements of discrete DMP
cs_r = CanonicalSystem( 'rhythmic', data_r.tau, 1.0 );
fs_r = NonlinearForcingTerm( cs_r, N );

% First, generate the trajectories
tscl_arr = [ 0.5, 1.0, 2.0 ];

% The Parameters for Forward Simulation
t0i   = 0.0;
T     = data_r.tau;
dt    = 1e-3;
t_arr = 0:dt:(T*2*pi*2);

y_data_arr = zeros( 2, length( t_arr ), length( tscl_arr ) );

scl = 1.8;
for i = 1 : length( tscl_arr )
    tscl = tscl_arr( i );

    % Rollout for data
    cs_tmp        = CanonicalSystem( 'rhythmic', tscl*data_r.tau, 1.0 );
    fs_tmp        = NonlinearForcingTerm( cs_tmp, N );
    trans_sys_tmp = TransformationSystem( data_r.alpha_z, data_r.beta_z, cs_tmp );
        
    % Calculate the nonlinear forcing term for discrete movement and rollout
    input_arr = fs_tmp.calc_forcing_term( t_arr( 1:end-1 ), data_r.weight, t0i, eye( 2 ) );
  
    [ y_arr, ~, ~ ] = trans_sys_tmp.rollout( zeros( 2, 1), zeros( 2, 1 ), zeros( 2, 1 ), scl*input_arr + data_r.alpha_z*data_r.beta_z*data_r.goal, t0i, t_arr  );  

    y_data_arr( :, :, i ) = y_arr;
    
end

ttmp_arr = data_r.tau * (0.25:0.25:2.0)*2*pi;
offset  = 66*( 1:length( ttmp_arr ) );

% Going through the iterations
for i = 1: length( tscl_arr )
    a = subplot( length( tscl_arr ), 2, 2*i, 'parent', f );
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
    title( a, ['$\tau=', num2str( tscl_arr( i ), '%.2f' ), '\tau^{(d)}$'], 'fontsize', 40 )
    
end

fig_save( f, './images/fig1c' )

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
