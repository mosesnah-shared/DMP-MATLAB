% [Title]     Contraction Theory with Position, Gradation of Movement
% [Author]    Moses Chong-ook Nah
% [Email]     mosesnah@mit.edu
% [Update]    At 2023.11.20

%% (--) Initialization
clear; close all; clc;
fig_config( 'fontSize', 20, 'markerSize', 10 )

%% =======================================================
%% (1-) Data Collection and Trajectory Generation
%%  -- (1A) Calling the trajectory

% The string list for the types of trajectories
% We have four trajectories for this example
traj_names = [ "cosine", "radial", "drawM_pos", "lift_up_down_pos" ];
Ntraj = length( traj_names );

traj_data = cell( 1, Ntraj );

% Loading the data
for i = 1 : Ntraj
    tmp = load( [ './learned_parameters/', traj_names{ i }, '.mat' ] );
    traj_data{ i } = tmp.data;
end

% Set the same duration, although will be generalized to different cases.
traj_data{ 3 }.tau = 3;
traj_data{ 4 }.tau = 3;

% Color array for plots
c_arr = [ 0.0000, 0.4470, 0.7410;
          0.8500, 0.3250, 0.0980;
          0.9290, 0.6940, 0.1250;
          0.4940, 0.1840, 0.5560;
          0.4660, 0.6740, 0.1880;
          0.3010, 0.7450, 0.9330];

clear tmp*

%%  -- (1B) Generating the trajectories 

t0i   =    0.0;           % The initial time of the movement rollout
T     =    5.0;           % The   whole time of the simulation 
dt    =   1e-3;           % Time-step  for the simulation
t_arr = 0:dt:T;           % Time array for the simulation

% Position and velocity array of the trajectory
 p_data_arr = zeros( 3, length( t_arr ), Ntraj );
dp_data_arr = zeros( 3, length( t_arr ), Ntraj );

% Scaling of the trajectory
scl_arr = [ 0.7, 1.0, 2.5, 5.0 ];

% The input force and goal arrays for the movements
force_arr = zeros( 3, length( t_arr )-1, 4 );
goal_arr  = zeros( 3, 4 );


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
    input_arr = fs.calc_forcing_term( t_arr( 1:end-1 ), data.weight, t0i, eye( 3 ) );

    % Scaling of Trajectory
    scl = scl_arr( i );
        
    % Calculate the new goal location and save
    new_goal = scl*data.goal;
    goal_arr( :, i ) = new_goal;

    % Calculate the new force and save
    new_force = scl*input_arr; 
    force_arr( :, :, i ) = new_force;

    % Rollout, note that initial condition is set to be zeros. 
    [ y_arr, z_arr, dy_arr ] = trans_sys.rollout( zeros( 3, 1 ), zeros( 3, 1 ), new_goal, new_force, t0i, t_arr  );    
    
     p_data_arr( :, :, i ) =  y_arr;
    dp_data_arr( :, :, i ) = dy_arr;

end

% Plotting the positions to double cehck
f = figure( ); a = axes( 'parent', f );
axis equal; view( 3 );
hold( a, 'on' )

for i = 1 : Ntraj
    plot3( a, p_data_arr( 1, :, i ), p_data_arr( 2, :, i ), p_data_arr( 3, :, i ), 'linewidth', 5, 'color', c_arr( i, : ) )
    scatter3( a, p_data_arr( 1,   1, i ), p_data_arr( 2,   1, i ), p_data_arr( 3,   1, i ), 'filled', 'o', 'markerfacecolor', 'white', 'markeredgecolor', c_arr( i, : ), 'linewidth', 4 )
    scatter3( a, p_data_arr( 1, end, i ), p_data_arr( 2, end, i ), p_data_arr( 3, end, i ), 'filled', 'd', 'markerfacecolor', 'white', 'markeredgecolor', c_arr( i, : ), 'linewidth', 4 )    
end

%%  -- (1C) Gradation of Movements

% Select the two movements that we will add
idx1 = 1;
idx2 = 4;

assert( any( idx1 == 1:4 ) && any( idx2 == 1:4 ) );

% The weight combination for both movements 
weight = 0:0.1:1.0;
Nw     = length( weight );

% Produce the position for these trajectory 
p_combined_arr = zeros( 3, length( t_arr ), Nw );

for i = 1 : Nw
    w = weight( i );

    % Get the force array
    f1_arr = force_arr( :, :, idx1 );
    f2_arr = force_arr( :, :, idx2 );

    % The goal array
    g1 = goal_arr( :, idx1 );
    g2 = goal_arr( :, idx2 );

    gnew =     g1 * w +     g2 * (1-w);
    fnew = f1_arr * w + f2_arr * (1-w);

    % The Three elements of DMP
    cs        = CanonicalSystem( 'discrete', data.tau, data.alpha_s );
    fs        = NonlinearForcingTerm( cs, N );
    trans_sys = TransformationSystem( data.alpha_z, data.beta_z, cs );

    % Rollout, note that initial condition is set to be zeros. 
    [ y_arr, ~, ~] = trans_sys.rollout( zeros( 3, 1 ), zeros( 3, 1 ), gnew, fnew, t0i, t_arr  );    
    
     p_combined_arr( :, :, i ) =  y_arr;

end

% Plot the movements

% Plotting the positions to double cehck
f = figure( ); a = axes( 'parent', f );
axis equal; view( 3 );
hold( a, 'on' )

c_arr2 = flipud( colormap( copper( 2*Nw ) ) );

for i = 1 : Nw
    plot3( a, p_combined_arr( 1, :, i ), p_combined_arr( 2, :, i ), p_combined_arr( 3, :, i ), 'linewidth', 5, 'color', c_arr2( i, : ) )
    scatter3( a, p_combined_arr( 1,   1, i ), p_combined_arr( 2,   1, i ), p_combined_arr( 3,   1, i ), 'filled', 'o', 'markerfacecolor', 'white', 'markeredgecolor', c_arr2( i, : ), 'linewidth', 4 )
    scatter3( a, p_combined_arr( 1, end, i ), p_combined_arr( 2, end, i ), p_combined_arr( 3, end, i ), 'filled', 'd', 'markerfacecolor', 'white', 'markeredgecolor', c_arr2( i, : ), 'linewidth', 4 )    
end
