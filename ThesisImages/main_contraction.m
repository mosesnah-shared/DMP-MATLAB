% [Title]     Generating Images for Contraction Theory Paper
% [Author]    Moses Chong-ook Nah
% [Email]     mosesnah@mit.edu
% [Update]    At 2024.02.12

%% (--) Initialization
clear; close all; clc;
fig_config( 'fontSize', 20, 'markerSize', 10 )

%% (1A) Dynamic Movement Primitives for multiple movements
close all
% For discrete movement, load the data
tmp  = load( './learned_parameters/discrete/A.mat' );
data = tmp.data;

% Getting the number of basis functions from the nonlinear forcing term
[ ~, N ] = size( data.weight );

cs1        = CanonicalSystem( 'discrete', data.tau, data.alpha_s );
fs1        = NonlinearForcingTerm( cs1, N );
trans_sys1 = TransformationSystem( data.alpha_z, data.beta_z, cs1 );

% The parameters for forward simulation
t0i   = 0.0;
T     = 16;
dt    = 1e-2;
t_arr = 0:dt:T;
Nt    = length( t_arr );

% Calculate the nonlinear forcing term for discrete movement
input_arr = fs1.calc_forcing_term( t_arr( 1:end-1 ), data.weight, t0i, eye( 2 ), 'trimmed' );

% Rollout, note that initial condition is set to be zeros. 
[ y_arr, z_arr, dy_arr ] = trans_sys1.rollout( zeros( 2, 1 ), zeros( 2, 1 ), data.goal, input_arr, t0i, t_arr  );  

f = figure( ); a = axes( 'parent', f );
a1 = subplot( 3, 2, 1 );
plot( a1, t_arr, cs1.calc( t_arr ), 'linewidth', 4, 'color',  [0.8500 0.3250 0.0980] );
set( a1, 'xticklabel', {}, 'yticklabel', {}, 'xlim', [0,16] )
xlabel( '$t$ (s)', 'fontsize', 40 )
ylabel( '$s_d$', 'fontsize', 40 )

a3 = subplot( 3, 2, 3 );
hold on
plot( a3, t_arr( 1:end-1 ), input_arr( 1, :), 'linewidth', 4, 'linestyle', '-' , 'color', [0.8500 0.3250 0.0980]	);
plot( a3, t_arr( 1:end-1 ), input_arr( 2, :), 'linewidth', 4, 'linestyle', '-' , 'color', [0.8500 0.3250 0.0980]	);
set( a3, 'xticklabel', {}, 'yticklabel', {}, 'xlim', [0,6.0], 'ylim', 0.3*[-1e+6,1e+6] )
xlabel( '$t$ (s)', 'fontsize', 40 )
ylabel( '$\mathbf{F}_d$', 'fontsize', 40 )

a5 = subplot( 3, 2, 5 );
hold on
plot( a5, y_arr( 1, : ), y_arr( 2, : ), 'linewidth', 4, 'color', [0.8500 0.3250 0.0980] );
scatter( a5, y_arr( 1,   1 ), y_arr( 2,   1 ), 500, 'filled',  'o', 'markerfacecolor', 'w', 'markeredgecolor', [0.8500 0.3250 0.0980], 'linewidth', 4 )
scatter( a5, y_arr( 1, end ), y_arr( 2, end ), 500, 'filled',  'd', 'markerfacecolor', 'w', 'markeredgecolor', [0.8500 0.3250 0.0980], 'linewidth', 4 )
set( a5, 'xticklabel', {}, 'yticklabel', {}, 'xlim', [-3, 9.5], 'ylim', [-2, 10.5] )
axis equal
xlabel( '$X$ (m)', 'fontsize', 40 )
ylabel( '$Y$ (m)', 'fontsize', 40 )

% For rhythmic movement, load the data
tmp  = load( './learned_parameters/rhythmic/heart.mat' );
data = tmp.data;

% Getting the number of basis functions from the nonlinear forcing term
[ ~, N ] = size( data.weight );

cs2        = CanonicalSystem( 'rhythmic', data.tau, 1.0 );
fs2        = NonlinearForcingTerm( cs2, N );
trans_sys2 = TransformationSystem( data.alpha_z, data.beta_z, cs2 );

% The parameters for forward simulation
t0i   = 0.0;
T     = 3;
dt    = 1e-4;
t_arr = 0:dt:T;
Nt    = length( t_arr );

input_arr2 = fs2.calc_forcing_term( t_arr( 1:end-1 ), data.weight , t0i, eye( 2 ) );
[ y_arr2, ~, ~ ] = trans_sys2.rollout( data.p_init, zeros( 2, 1 ), [ 0; 0.8 ], input_arr2, t0i, t_arr  );

a2 = subplot( 3, 2, 2 );
plot( t_arr, cs2.calc( t_arr ), 'linewidth', 4, 'color', [0 0.4470 0.7410] );
set( a2, 'xticklabel', {}, 'yticklabel', {}, 'xlim', [0,T] )
xlabel( '$t$ (s)', 'fontsize', 40 )
ylabel( '$s_r$', 'fontsize', 40 )

a4 = subplot( 3, 2, 4 );
hold on
plot( a4, t_arr( 1:end-1 ), input_arr2( 1, :), 'linewidth', 4, 'linestyle', '-' , 'color', [0 0.4470 0.7410] );
plot( a4, t_arr( 1:end-1 ), input_arr2( 2, :), 'linewidth', 4, 'linestyle', '-' , 'color', [0 0.4470 0.7410] );
set( a4, 'xticklabel', {}, 'yticklabel', {}, 'xlim', [0,T], 'ylim', [-6*10^5,6*10^5] )
xlabel( '$t$ (s)', 'fontsize', 40 )
ylabel( '$\mathbf{F}_r$', 'fontsize', 40 )

a6 = subplot( 3, 2, 6 );

plot( a6, y_arr2( 1, 1:end-1 ), y_arr2( 2, 1:end-1 ), 'linewidth', 4, 'color', [0 0.4470 0.7410] );
set( a6, 'xticklabel', {}, 'yticklabel', {}, 'xlim', [-20,20], 'ylim',[-20,20])
xlabel( '$X$ (m)', 'fontsize', 40 )
ylabel( '$Y$ (m)', 'fontsize', 40 )
axis equal

export_fig ./Images/DMP_images/dis_and_rhy_DMP.pdf -transparent

%% (1B) Dynamic Movement Primitives for Sequencing Movements, Data
%% (--) (-1a) Data Processing
close all
% Import the learned weights
alphabets = { 'A_loose', 'B_loose', 'C_loose' };
Na = length( alphabets );

% dataset 
traj_data = cell( 1, Na );

for i = 1 : Na
    a = alphabets{ i };

    % Call Data
    load( [ './learned_parameters/discrete/', a, '.mat' ] );
    traj_data{ i } = data;
end

Ntraj = Na;

tinit =    1.0;           % The initial time of the simulation
T     =   30.0;           % The   whole time of the simulation 
toff  =   -0.1;           % Time-offset for the subsequent movement to start, must be negative
dt    =   1e-3;           % Time-step  for the simulation
t_arr = 0:dt:T;           % Time array for the simulation

% Time offset must be strictly negative. 
assert( toff < 0 );

% Position and velocity array of the trajectory
% This is 2D Data
 p_data_arr = zeros( 2, length( t_arr ), Ntraj );
dp_data_arr = zeros( 2, length( t_arr ), Ntraj );

% Scaling and Rotation of the trajectory
scl_arr = [ 1.3, 1.2, 1.1 ];

% The input force array for the movements
force_arr = zeros( 2, length( t_arr )-1, Ntraj );

% The goal and initated time and final time for each movement
goal_arr = zeros( 2, Ntraj );
 t0i_arr = zeros( 1, Ntraj );
 t0f_arr = zeros( 1, Ntraj );

% The first movement starts at tinit, hence setting up
 t0i = tinit;          

% We need to first generate the DMP 
for i = 1 : Ntraj
    
    % Call the data from the trajectory
    data = traj_data{ i };

    % Saving the initial time in the array
    t0i_arr( i ) = t0i;
    t0f_arr( i ) = t0i + data.tau;
    
    % Getting the number of basis functions from the nonlinear forcing term
    [ ~, N ] = size( data.weight );

    % The Three elements of DMP
    cs        = CanonicalSystem( 'discrete', data.tau, data.alpha_s );
    fs        = NonlinearForcingTerm( cs, N );
    trans_sys = TransformationSystem( data.alpha_z, data.beta_z, cs );

    % Calculate the nonlinear forcing term 
    % This can be diminishing
    input_arr_trimmed = fs.calc_forcing_term( t_arr( 1:end-1 ), data.weight, t0i, eye( 2 ), 'trimmed' );
    input_arr         = fs.calc_forcing_term( t_arr( 1:end-1 ), data.weight, t0i, eye( 2 ) );

    % Scaling
    scl = scl_arr( i );
        
    % Calculate the new goal location and save
    new_goal = scl*data.goal;
    goal_arr( :, i ) = new_goal;

    % Calculate the new force and save
    new_force = scl*input_arr; 
    new_force_trimmed = scl*input_arr_trimmed; 

    if i == Ntraj
        force_arr( :, :, i ) = new_force_trimmed;
    else
        force_arr( :, :, i ) = new_force;
    end

    % Rollout, note that initial condition is set to be zeros. 
    [ y_arr, z_arr, dy_arr ] = trans_sys.rollout( zeros( 2, 1 ), zeros( 2, 1 ), new_goal, new_force_trimmed, t0i, t_arr  );    
    
     p_data_arr( :, :, i ) =  y_arr;
    dp_data_arr( :, :, i ) = dy_arr;

    % Define the next movement's starting time.
    t0i = t0i + toff + data.tau;

end

% Add position offset
xy_off = [ 0,0 ; 11.5, 10.5; 30, 9 ]';

p_data_arr( :, :, 2 ) = p_data_arr( :, :, 2 ) + xy_off( :, 2 );
p_data_arr( :, :, 3 ) = p_data_arr( :, :, 3 ) + xy_off( :, 3 );

goal_arr( :, 1 ) = p_data_arr( :, end, 1 );
goal_arr( :, 2 ) = p_data_arr( :, end, 2 );
goal_arr( :, 3 ) = p_data_arr( :, end, 3 );

% Plotting the positions to double cehck
f = figure( ); a = axes( 'parent', f );
axis equal; 
hold( a, 'on' )

for i = 1 : Ntraj
    plot( a, p_data_arr( 1, :, i ), p_data_arr( 2, :, i ), 'linewidth', 5, 'color', 'k' )
    scatter( a, p_data_arr( 1,   1, i ), p_data_arr( 2,   1, i ), 300, 'filled', 'o', 'markerfacecolor', 'white', 'markeredgecolor', 'k', 'linewidth', 4 )
    scatter( a, p_data_arr( 1, end, i ), p_data_arr( 2, end, i ), 300, 'filled', 'd', 'markerfacecolor', 'white', 'markeredgecolor', 'k', 'linewidth', 4 )    
end

% Setting the Departure and Arrival Time 
t_depart_arr = t0f_arr( 1:2 ) - [ 0.4, 0.1 ];
t_arrive_arr = t0i_arr( 2:3 ) + [ 0.0, 0.2 ];

% Define the A matrix and az, bz values
% All should be identical, 
% Tau can be arbitrarily chosen based on the proof
% Choosing one arbitrarily
tau = traj_data{ 1 }.tau;
as  = traj_data{ 1 }.alpha_s;
az  = traj_data{ 1 }.alpha_z;
bz  = traj_data{ 1 }.beta_z;

% A Matrix
Amat = 1/tau * [            -as,       zeros( 1, 2 ),  zeros( 1, 2 ); 
                  zeros( 2, 1 ),       zeros( 2, 2 ),       eye( 2 );
                  zeros( 2, 1 ), -az * bz * eye( 2 ), -az * eye( 2 )];

Cg = 1/tau * [ zeros( 3, 1 ); az * bz * goal_arr( :, 1 )];

% The initial condition is identical to the first trajectory x0 value.
x0 = [ 1; p_data_arr( :, 1, 1 ); dp_data_arr( :, 1, 1 )*traj_data{1}.tau];

% The coupled x array 
xc_arr = zeros( 5, length( t_arr ) );
xc_arr( :, 1 ) = x0;

xcurr = x0;
% Iterating through the x_arr

% Integrating over time. 
for i = 1 : length( t_arr )-1
    
    t = t_arr( i );

    if t <= tinit
        dx = zeros( 5, 1 );
    else
        % Add the force for the first movement
        dx  = Amat * xcurr + 1/traj_data{ 1 }.tau * [ zeros( 3, 1 ); force_arr( :, i, 1 ) ] + Cg;

        % Iterating through the movement with activation
        for j = 2 : Ntraj
            fp = force_arr( :, i, j-1 );  % Force for previous movement
            fn = force_arr( :, i, j   );  % Force for the next movement
            
            Dp = traj_data{ j-1 }.tau;
            Dn = traj_data{ j   }.tau;

            td   = t_depart_arr( j-1 );
            ta   = t_arrive_arr( j-1 );
        
            gp = 1/Dp * [ zeros( 3, 1 ); az * bz * goal_arr( :, j-1 ) ];
            gn = 1/Dn * [ zeros( 3, 1 ); az * bz * goal_arr( :, j   ) ];
            
            gain = clip_func( t, td, ta );
            dx  = dx +  gain * ( - 1/Dp * [ zeros( 3, 1 ); fp ] + 1/Dn * [ zeros( 3, 1 ); fn ] ) ...
                     +  gain * ( -gp + gn ); 
        end

    end

    xc_arr( :, i+1 ) = xc_arr( :, i ) + dx * dt;     
    xcurr = xc_arr( :, i+1 );

end

% plot( a, xc_arr( 2, : ), xc_arr( 3, : ))
close all
% Drawing the Figure 
% First, the first letter 
f = figure( ); a = axes( 'parent', f );
hold on

cblue   = [0 0.4470 0.7410];
corange = [0 0.4470 0.7410];

set( a , 'xlim', [-5, 35 ], 'ylim', [-6, 16], 'xticklabel', {}, 'yticklabel', {})
plot( a, p_data_arr( 1, :, 1 ), p_data_arr( 2, :, 1 ), 'linewidth', 5, 'color', 0.7*ones( 1, 3 ) )

plot( a, p_data_arr( 1, :, 2 ), p_data_arr( 2, :, 2 ), 'linewidth', 5, 'color', 0.7*ones( 1, 3 ) )
scatter( a, p_data_arr( 1,   1, 2 ), p_data_arr( 2,   1, 2 ), 300, 'filled', 'o', 'markerfacecolor', 'white', 'markeredgecolor', 'k', 'linewidth', 4 )
scatter( a, p_data_arr( 1, end, 2 ), p_data_arr( 2, end, 2 ), 300, 'filled', 'd', 'markerfacecolor', 'white', 'markeredgecolor', 'k', 'linewidth', 4 )    

plot( a, p_data_arr( 1, :, 3 ), p_data_arr( 2, :, 3 ), 'linewidth', 5, 'color', 0.7*ones( 1, 3 ) )
scatter( a, p_data_arr( 1,   1, 3 ), p_data_arr( 2,   1, 3 ), 300, 'filled', 'o', 'markerfacecolor', 'white', 'markeredgecolor', 'k', 'linewidth', 4 )
scatter( a, p_data_arr( 1, end, 3 ), p_data_arr( 2, end, 3 ), 300, 'filled', 'd', 'markerfacecolor', 'white', 'markeredgecolor', 'k', 'linewidth', 4 )    

N1 = 5000;
plot( a, xc_arr( 2, 1:N1 ), xc_arr( 3, 1:N1 ), 'linewidth', 7, 'color', cblue )
scatter( a, xc_arr( 2, N1 ), xc_arr( 3, N1 ), 600, 'filled', 'o', 'markerfacecolor', 'white', 'markeredgecolor', cblue, 'linewidth', 4 )
scatter( a, p_data_arr( 1,   1, 1 ), p_data_arr( 2,   1, 1 ), 300, 'filled', 'o', 'markerfacecolor', 'white', 'markeredgecolor', cblue, 'linewidth', 4 )
scatter( a, p_data_arr( 1, end, 1 ), p_data_arr( 2, end, 1 ), 300, 'filled', 'd', 'markerfacecolor', 'white', 'markeredgecolor', 'k', 'linewidth', 4 )    

export_fig ./Images/DMP_images/sequenceDMP_1.pdf -transparent

%% (--) (-1c) Plotting, Figure 2
close all
f = figure( ); a = axes( 'parent', f );
hold on

cblue   = [0 0.4470 0.7410];
corange = [0 0.4470 0.7410];

set( a , 'xlim', [-5, 35 ], 'ylim', [-6, 16], 'xticklabel', {}, 'yticklabel', {})
plot( a, p_data_arr( 1, :, 1 ), p_data_arr( 2, :, 1 ), 'linewidth', 5, 'color', cblue )

plot( a, p_data_arr( 1, :, 2 ), p_data_arr( 2, :, 2 ), 'linewidth', 5, 'color', 0.7*ones( 1, 3 ) )
scatter( a, p_data_arr( 1,   1, 2 ), p_data_arr( 2,   1, 2 ), 300, 'filled', 'o', 'markerfacecolor', 'white', 'markeredgecolor', 'k', 'linewidth', 4 )
scatter( a, p_data_arr( 1, end, 2 ), p_data_arr( 2, end, 2 ), 300, 'filled', 'd', 'markerfacecolor', 'white', 'markeredgecolor', cblue, 'linewidth', 4 )    

plot( a, p_data_arr( 1, :, 3 ), p_data_arr( 2, :, 3 ), 'linewidth', 5, 'color', 0.7*ones( 1, 3 ) )
scatter( a, p_data_arr( 1,   1, 3 ), p_data_arr( 2,   1, 3 ), 300, 'filled', 'o', 'markerfacecolor', 'white', 'markeredgecolor', 'k', 'linewidth', 4 )
scatter( a, p_data_arr( 1, end, 3 ), p_data_arr( 2, end, 3 ), 300, 'filled', 'd', 'markerfacecolor', 'white', 'markeredgecolor', 'k', 'linewidth', 4 )    

N1 = 12500;
plot( a, xc_arr( 2, 1:N1 ), xc_arr( 3, 1:N1 ), 'linewidth', 7, 'color', cblue )
scatter( a, xc_arr( 2, N1 ), xc_arr( 3, N1 ), 600, 'filled', 'o', 'markerfacecolor', 'white', 'markeredgecolor', cblue, 'linewidth', 4 )
scatter( a, p_data_arr( 1,   1, 1 ), p_data_arr( 2,   1, 1 ), 300, 'filled', 'o', 'markerfacecolor', 'white', 'markeredgecolor', cblue, 'linewidth', 4 )
scatter( a, p_data_arr( 1, end, 1 ), p_data_arr( 2, end, 1 ), 300, 'filled', 'd', 'markerfacecolor', 'white', 'markeredgecolor', cblue, 'linewidth', 4 )    
export_fig ./Images/DMP_images/sequenceDMP_2.pdf -transparent

%% (--) (-1d) Plotting, Figure 3
close all
f = figure( ); a = axes( 'parent', f );
hold on

cblue   = [0 0.4470 0.7410];
corange = [0 0.4470 0.7410];

set( a , 'xlim', [-5, 35 ], 'ylim', [-6, 16], 'xticklabel', {}, 'yticklabel', {})
plot( a, p_data_arr( 1, :, 1 ), p_data_arr( 2, :, 1 ), 'linewidth', 5, 'color', cblue )

plot( a, p_data_arr( 1, :, 2 ), p_data_arr( 2, :, 2 ), 'linewidth', 5, 'color', cblue )
scatter( a, p_data_arr( 1,   1, 2 ), p_data_arr( 2,   1, 2 ), 300, 'filled', 'o', 'markerfacecolor', 'white', 'markeredgecolor', cblue, 'linewidth', 4 )
scatter( a, p_data_arr( 1, end, 2 ), p_data_arr( 2, end, 2 ), 300, 'filled', 'd', 'markerfacecolor', 'white', 'markeredgecolor', cblue, 'linewidth', 4 )    

plot( a, p_data_arr( 1, :, 3 ), p_data_arr( 2, :, 3 ), 'linewidth', 5, 'color', 0.7*ones( 1, 3 ) )
scatter( a, p_data_arr( 1,   1, 3 ), p_data_arr( 2,   1, 3 ), 300, 'filled', 'o', 'markerfacecolor', 'white', 'markeredgecolor', cblue, 'linewidth', 4 )
scatter( a, p_data_arr( 1, end, 3 ), p_data_arr( 2, end, 3 ), 300, 'filled', 'd', 'markerfacecolor', 'white', 'markeredgecolor', 'k', 'linewidth', 4 )    

N1 = 20500;
plot( a, xc_arr( 2, 1:N1 ), xc_arr( 3, 1:N1 ), 'linewidth', 7, 'color', cblue )
scatter( a, xc_arr( 2, N1 ), xc_arr( 3, N1 ), 600, 'filled', 'o', 'markerfacecolor', 'white', 'markeredgecolor', cblue, 'linewidth', 4 )
scatter( a, p_data_arr( 1,   1, 1 ), p_data_arr( 2,   1, 1 ), 300, 'filled', 'o', 'markerfacecolor', 'white', 'markeredgecolor', cblue, 'linewidth', 4 )
scatter( a, p_data_arr( 1, end, 1 ), p_data_arr( 2, end, 1 ), 300, 'filled', 'd', 'markerfacecolor', 'white', 'markeredgecolor', cblue, 'linewidth', 4 )    

export_fig ./Images/DMP_images/sequenceDMP_3.pdf -transparent


%% (1C) Dynamic Movement Primitives for Synchronized Movements
%% (--) (-1a) Data Processing

% Doing it for a 6DOF Robot
Nosc = 6;

c_arr = cell( 1, 7 );
c_arr{ 1 } = [0.0000 0.4470 0.7410];
c_arr{ 2 } = [0.8500 0.3250 0.0980];
c_arr{ 3 } = [0.9290 0.6940 0.1250];
c_arr{ 4 } = [0.4940 0.1840 0.5560];
c_arr{ 5 } = [0.4660 0.6740 0.1880];	
c_arr{ 6 } = [0.3010 0.7450 0.9330];
c_arr{ 7 } = [0.6350 0.0780 0.1840];

% Now, define three complex numbers for the three oscillators 
% Parameters of the Andronov-Hopf Oscillator 
r = 1.0; w = 2*pi;

% Using complex number
r0 = 0.5;

r_arr = 0.9*ones( 1, Nosc );
phase_off = 2*pi*rand( 1, Nosc )-2*pi;

% First, run the simulation for both oscillators
dt = 1e-4;
T  = 6.0;
N  = round( T/dt );

% Time array 
t_arr = (0:N)*dt;
z_arr = zeros( Nosc, N+1 );
z_arr( :, 1 ) = r_arr .* exp( 1i * ( pi/2 + phase_off ) );

% This is a zero matrix without any synchronization
M1 = zeros( Nosc, Nosc );

k = 6.0;
M2 = diag( -k*ones( 1, Nosc ) ) +k * circshift( diag( exp( 1i* (2*pi/Nosc) ) * ones( 1, Nosc )  ) ,1,2  );

M = M1;

for i = 2 : N+1

    % Calculating dz 
    dz = zeros( 6, 1 );
    z_old = z_arr( :, i-1 );
    for j = 1 : Nosc
        dz( j ) =  ( ( r^2-abs( z_old( j ) )^2 ) + 1i*w ) * z_old( j ) + M( j, : ) * z_old;
    end

    if i == round( N/8 )
        M=M2;
    end

    z_arr( :, i ) = z_old + dz*dt;
end

%% (--) (-1ab) Drawing the Main Plots

f = figure( ); a = axes( 'parent' ,f );
hold on

for i = 1 : Nosc
    plot( a, t_arr, angle( z_arr( i, : ) ), 'linewidth', 5, 'color', c_arr{ i } );
end
set( a, 'xticklabel', {}, 'yticklabel', {} )
xlabel( a, '$t$ (s)', 'fontsize', 40)

export_fig ./Images/synchronization/sync1.pdf -transparent


%% (--) (-1b) Drawing the Plots
close all;

f = figure( ); 

%Setting the moment
Nm = 6000;
for i = 1 : Nosc
    a = subplot( 1, Nosc, i  ); set( a, 'parent', f );
    hold on;
    plot( a, r*cos( 0:0.001:2*pi),r*sin( 0:0.001:2*pi), 'linewidth', 5, 'color', c_arr{ i } )
    scatter( a, real( z_arr( i, Nm ) ), imag( z_arr( i, Nm ) ), 500, 'o', 'filled', 'linewidth', 5, 'markeredgecolor', c_arr{ i }, 'markerfacecolor', 'w' );
    axis equal
    set( a, 'xlim', [-2, 2], 'ylim', [-2, 2], 'xticklabel', {}, 'yticklabel', {})
end


export_fig ./Images/synchronization/sync2.pdf -transparent

f = figure( ); 
%Setting the moment
Nm = 50000;
for i = 1 : Nosc
    a = subplot( 1, Nosc, i  ); set( a, 'parent', f );
    hold on;
    plot( a, r*cos( 0:0.001:2*pi),r*sin( 0:0.001:2*pi), 'linewidth', 5, 'color', c_arr{ i } )
    scatter( a, real( z_arr( i, Nm ) ), imag( z_arr( i, Nm ) ), 500, 'o', 'filled', 'linewidth', 5, 'markeredgecolor', c_arr{ i }, 'markerfacecolor', 'w' );
    axis equal
    set( a, 'xlim', [-2, 2], 'ylim', [-2, 2], 'xticklabel', {}, 'yticklabel', {})
end

export_fig ./Images/synchronization/sync3.pdf -transparent

%% (--) (-1c) Drawing the Plots
% With this, generating a robot motion from this 
% Getting the angle of z_arr 

ang = angle( z_arr );

% We want to generate a motion of 
% ang, ang +2pi*1/6, ang+2pi*2/6,...

% Using cosine for the angles 
joint_angs = cos( ang );

% Get x and y values
x = cos( joint_angs + pi/2 );
y = sin( joint_angs + pi/2 );

xorig = zeros( 1, length( t_arr ) );
yorig = zeros( 1, length( t_arr ) );
x_abs = cumsum( x, 1 );
y_abs = cumsum( y, 1 );

x_abs = [ xorig; x_abs ];
y_abs = [ yorig; y_abs ];

% Moments
N_arr = round( N/8) - 600*(1:11)-500;
alpha = 0.5*ones( 1, length(N_arr));
f = figure(); a = axes( 'parent', f );
hold on

for i = 1 : Nosc
    for j = 1 : length( N_arr )
        idx = N_arr( j );
        xtmp = x_abs( i:i+1, idx );
        ytmp = y_abs( i:i+1, idx );
        plot( a, xtmp, ytmp, 'linewidth', 5, 'color', [c_arr{ i }, alpha( j )] );

        xtmp = x_abs( :, idx );
        ytmp = y_abs( :, idx);        
        scatter( a, xtmp( i ), ytmp( i ), 500, 'o', 'filled', 'markerfacecolor', 'w', 'markeredgecolor', c_arr{ i }, 'linewidth', 5 )        

        if i == Nosc
            % Final end-effector
            scatter( a, xtmp( end ), ytmp( end ), 500, 'o', 'filled', 'markerfacecolor', 'w', ...
                                'markeredgecolor', 'k', 'markeredgealpha', alpha( j ), 'markerfacealpha', alpha( j ), 'linewidth', 5 )
        end        
    end

end

for i = 1 : Nosc
    plot( a, x_abs( i, min( N_arr ):max( N_arr ) ), y_abs( i, min( N_arr ):max( N_arr ) ), 'linewidth', 5, 'color', c_arr{ i } );
end
plot( a, x_abs( end, min( N_arr ):max( N_arr ) ), y_abs( end, min( N_arr ):max( N_arr ) ), 'linewidth', 5, 'color', 'k' );

axis equal
set( a, 'xlim', [-3.5, 3.5], 'ylim', [-1, 6], 'xticklabel', {}, 'yticklabel', {})

export_fig ./Images/synchronization/sync4.pdf -transparent


%% (--) (-1b) Drawing the Plots
close all;

f = figure( ); 
for i = 1 : Nosc
    a = subplot( 1, Nosc, i  ); set( a, 'parent', f );
    hold on;
    plot( a, real( z_arr( i, : ) ), imag( z_arr( i, : ) ), 'linewidth', 5, 'color', c_arr{ i } );
    axis equal
    set( a, 'xlim', [-2, 2], 'ylim', [-2, 2])
end

% With this, generating a robot motion from this 
% Getting the angle of z_arr 

ang = angle( z_arr );

% We want to generate a motion of 
% ang, ang +2pi*1/6, ang+2pi*2/6,...

% Using cosine for the angles 
joint_angs = cos( ang );

% Get x and y values
x = cos( joint_angs + pi/2 );
y = sin( joint_angs + pi/2 );

xorig = zeros( 1, length( t_arr ) );
yorig = zeros( 1, length( t_arr ) );
x_abs = cumsum( x, 1 );
y_abs = cumsum( y, 1 );

x_abs = [ xorig; x_abs ];
y_abs = [ yorig; y_abs ];

% Moments
N_arr = round( N/8) + 1000*(1:11)+20000;
alpha = 0.5*ones( 1, length(N_arr));
f = figure(); a = axes( 'parent', f );
hold on

for i = 1 : Nosc
    for j = 1 : length( N_arr )
        idx = N_arr( j );
        xtmp = x_abs( i:i+1, idx );
        ytmp = y_abs( i:i+1, idx );
        pline = plot( a, xtmp, ytmp, 'linewidth', 5, 'color', [c_arr{ i }, alpha( j )] );

        xtmp = x_abs( :, idx );
        ytmp = y_abs( :, idx);        
        scatter( a, xtmp( i ), ytmp( i ), 500, 'o', 'filled', 'markerfacecolor', 'w', 'markeredgecolor', c_arr{ i }, 'linewidth', 5 )        

        if i == Nosc
            % Final end-effector
            scatter( a, xtmp( end ), ytmp( end ), 500, 'o', 'filled', 'markerfacecolor', 'w', ...
                                'markeredgecolor', 'k', 'markeredgealpha', alpha( j ), 'markerfacealpha', alpha( j ), 'linewidth', 5 )
        end        
    end


end


for i = 1 : Nosc
    plot( a, x_abs( i, min( N_arr ):max( N_arr ) ), y_abs( i, min( N_arr ):max( N_arr ) ), 'linewidth', 5, 'color', c_arr{ i } );
end
plot( a, x_abs( end, min( N_arr ):max( N_arr ) ), y_abs( end, min( N_arr ):max( N_arr ) ), 'linewidth', 5, 'color', 'k' );

axis equal
set( a, 'xlim', [-3.5, 3.5], 'ylim', [-1, 6], 'xticklabel', {}, 'yticklabel', {})

export_fig ./Images/synchronization/sync5.pdf -transparent


%% (1D) Dynamic Movement Primitives Synchronization Primitives


