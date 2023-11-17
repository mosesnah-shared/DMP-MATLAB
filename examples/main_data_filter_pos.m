% [Title]     Data Filter and learn the Weight Array, Position
% [Author]    Moses Chong-ook Nah
% [Email]     mosesnah@mit.edu
% [Update]    At 2023.11.13

%% (--) Initialization
clear; close all; clc;
fig_config( 'fontSize', 20, 'markerSize', 10 )

%% =======================================================
%% (1-) Position Data Filtering and Learning the Weights
%%  -- (1A) Calling the data and get the Forward Kinematics 

% All dataset are saved under 'data' directory
idx = 1;
traj_names = [ "lift_up_down", "drawM" ]; 
file_names = [ "example3/iiwa_example_orient", 'example4/iiwa_example_pos' ];


raw_data = parse_txt( [ 'data/', file_names{ idx }, '.txt' ], 0 );

% Read the time (with offset) and joint-trajectory data
t_arr_demo = raw_data( :, 1 )' - raw_data( 1, 1 );
q_arr_demo = raw_data( :, 2:8 )';

% Call the robot for Forward Kinematics
robot = iiwa14( 'high' );
robot.init( );

% The number of sample points P, which we denote as Np
[ ~, Np ] = size( q_arr_demo );

% End-effect Forward Kinematics, both position and orientation
% We will learn the position, but also saving the orientation for backup
p_arr_raw = zeros( 3, Np );
R_arr_raw = zeros( 3, 3, Np );

% Calculating the position/orientation with Forward Kinematics
for i = 1 : Np
    H = robot.getForwardKinematics( q_arr_demo( :, i ) );
    p_arr_raw( :, i )    = H( 1:3,   4 );
    R_arr_raw( :, :, i ) = H( 1:3, 1:3 ); 
end

% Minus the offset for initial dataset zero 
p_arr_raw = p_arr_raw - p_arr_raw( :, 1 );


f = figure( ); a = axes( 'parent', f );
hold( a, 'on' ); axis equal
plot3( a, p_arr_raw( 1, : ), p_arr_raw( 2, : ), p_arr_raw( 3, : ), 'linewidth', 3 )
scatter3( a, p_arr_raw( 1,   1 ), p_arr_raw( 2,   1 ), p_arr_raw( 3,   1 ), 500, 'filled', 'd', 'markerfacecolor', 'white', 'markeredgecolor', [0.0000 0.4470 0.7410], 'linewidth', 5)
scatter3( a, p_arr_raw( 1, end ), p_arr_raw( 2, end ), p_arr_raw( 3, end ), 500, 'filled', 'o', 'markerfacecolor', 'white', 'markeredgecolor', [0.8500 0.3250 0.0980]	, 'linewidth', 5)


%%  -- (1B) Filtering out the Position data

% Simple time diff for comparison
dp_arr_raw  = data_diff(  p_arr_raw );
ddp_arr_raw = data_diff( dp_arr_raw );

% Filtering out the Data
  p_arr_filt = zeros( size( p_arr_raw ), 'like', p_arr_raw );
 dp_arr_filt = zeros( size( p_arr_raw ), 'like', p_arr_raw );
ddp_arr_filt = zeros( size( p_arr_raw ), 'like', p_arr_raw );

% Filtered position, velocity and acceleration
for i = 1 : 3
    % Gaussian Filtering of Position 
     p_arr_filt( i, : ) = smoothdata( p_arr_raw( i, : ), "gaussian", 50 );

    % Diff and Gaussian Filter of Velocity  
   tmp_vel = data_diff( p_arr_filt );
    dp_arr_filt( i, : ) = smoothdata( tmp_vel( i, : ), "gaussian", 50 );

    % Diff and Gaussian Filter of Acceleration
   tmp_acc = data_diff( dp_arr_filt );
   ddp_arr_filt( i, : ) = smoothdata( tmp_acc( i, : ), "gaussian", 50 );   
end


% Plotting the Data
tmp_plot = cell( 1, 6 );
tmp_plot{ 1 } =   p_arr_raw; tmp_plot{ 2 } =   p_arr_filt;
tmp_plot{ 3 } =  dp_arr_raw; tmp_plot{ 4 } =  dp_arr_filt;
tmp_plot{ 5 } = ddp_arr_raw; tmp_plot{ 6 } = ddp_arr_filt;

tmp_label = [ "Pos.", "Vel.", "Acc." ];

% Row
for i = 1:3

    % Column
    for j = 1 :2
        subplot( 3, 2, 2*(i-1)+j )
        hold on
        plot( t_arr_demo, tmp_plot{ 2*(i-1)+j }, 'linewidth', 3 )
        set( gca, 'xlim', [ 0, max( t_arr_demo ) ], 'xticklabel', {}, 'fontsize', 30 )
        
        if j == 1
           ylabel( tmp_label( i ) )
           tmp = get( gca, 'ylim' );
        else
           set( gca, 'ylim', tmp ) 
        end
    end

end


%%  -- (1C) Learning the Weights for Imitation Learning

% Parameters of DMP
alpha_s   = 1.0;
alpha_z   = 1000.0;
beta_z    = 0.5 * alpha_z;
N         = 100;
tau       = max( t_arr_demo );

% The Three elements of DMP
cs        = CanonicalSystem( 'discrete', tau, alpha_s );
trans_sys = TransformationSystem( alpha_z, beta_z, cs );
fs        = NonlinearForcingTerm( cs, N );

% Calculating the required Nonlinear Forcing Term
P  = length( t_arr_demo );
y0d = p_arr_filt( :, 1   );     % Initial position of demonstration
gd  = p_arr_filt( :, end );     % Goal    position of demonstration
f_arr = trans_sys.get_desired( p_arr_filt, dp_arr_filt, ddp_arr_filt, gd );

% The phi matrix 
Phi_mat = zeros( P, N );

% Calculate the Phi matrix for Weight Learning
% Note that compared to Ijspeert 2013, this followed Koutras and Doulgeri (2020)
for i = 1 : P 
    t = t_arr_demo( i );
    Phi_mat( i, : ) = fs.calc_multiple_ith( t, 1:N )/ fs.calc_whole_at_t( t ) * cs.calc( t );
end

% Learning the weights with Linear Least-square fitting
w_arr = transpose( ( Phi_mat' * Phi_mat )^(-1) * Phi_mat' * f_arr' );


% Double check the Learned weights with basic DMP
% Initial conditions
y0 =  p_arr_filt( :,   1 );
z0 = dp_arr_filt( :,   1 )/tau; 
g  =  p_arr_filt( :, end );
t0i = 0.0;
T   = max( t_arr_demo );
t_arr = linspace( 0, T, 2000 );

% Calculate the nonlinear forcing term
input_arr = fs.calc_forcing_term( t_arr( 1:end-1 ), w_arr, t0i, eye( length( g ) ) );

[ y_arr, z_arr, dy_arr ] = trans_sys.rollout( y0, z0, g, input_arr, t0i, t_arr  );
    
f = figure( ); a = axes( 'parent', f );
hold( a, 'on' );
plot( a, t_arr_demo, p_arr_filt, 'linewidth', 10, 'linestyle', '--', 'color', 'black' )
plot( a,      t_arr,      y_arr, 'linewidth',  5, 'linestyle',  '-', 'color', 'blue'  )


%%  -- (1D) Generalize the demonstrated trajectory, under Rotation

% Since the weights are learned, we can scale and rotate the trajectory 
y0 = zeros( 3, 1 );
z0 = dp_arr_filt( :, 1 )/cs.tau;

% Parameters
t0i   = 0.5;
T     = 5.0;
t_arr = 0:1e-3:T;

f = figure( ); a = axes( 'parent', f );
hold( a, 'on' ); axis equal
set( a, 'view', [46.1836, 38.4642], 'visible', 'off' )
scatter3( a, 0, 0, 0, 200, 'filled', 'markeredgecolor', 'black',  'markerfacecolor', 'white', 'linewidth', 2 );

theta_arr = 1:60:360;

for i = 1 : length( theta_arr )

    theta = theta_arr( i ); 
    Rmat  = rotz( theta );
     
    % Goal Location and its rotation
    g  = p_arr_filt( :, end ) - p_arr_filt( :, 1 ); 
    g  = Rmat * g;
     
    input_arr = fs.calc_forcing_term( t_arr( 1:end-1), w_arr, t0i, eye( 3 ) );
    
    [ y_arr, z_arr, dy_arr ] = trans_sys.rollout( y0, z0, g, Rmat * input_arr, t0i, t_arr  );
    
    % The original demonstrated trajectory
    if i == 1 
        plot3( a, y_arr( 1, : ), y_arr( 2, : ), y_arr( 3, : ), 'linewidth', 5, 'color', [0.8500 0.3250 0.0980] )
    else
        plot3( a, y_arr( 1, : ), y_arr( 2, : ), y_arr( 3, : ), 'linewidth', 5, 'color', 'black' )
    end

    scatter3( a, g( 1 ), g( 2 ), g( 3 ), 200, 'filled', 'markeredgecolor', [0 0.4470 0.7410], 'markerfacecolor', 'white', 'linewidth', 4 );

end


%%  -- (1E) Generalize the demonstrated trajectory, under Scaling

% Since the weights are learned, we can scale and rotate the trajectory 
y0 = zeros( 3, 1 );
z0 = dp_arr_filt( :, 1 )/cs.tau;

% Parameters
t0i   = 0.0;
T     = 5.0;
t_arr = 0:1e-3:T;

f = figure( ); a = axes( 'parent', f );
hold( a, 'on' ); axis equal
set( a, 'view', [46.1836, 38.4642], 'visible', 'off' )
scatter3( a, 0, 0, 0, 200, 'filled', 'markeredgecolor', 'black',  'markerfacecolor', 'white', 'linewidth', 2 );

scl_arr = [ 1.0, 1.2, 1.5, 2.0, 4.0];

for i = 1 : length( scl_arr )

    scl = scl_arr( i ); 
     
    % Goal Location and its rotation
    g  = p_arr_filt( :, end ) - p_arr_filt( :, 1 ); 
    g  = scl * g;
     
    input_arr = fs.calc_forcing_term( t_arr( 1:end-1), w_arr, t0i, eye( 3 ) );
    
    [ y_arr, z_arr, dy_arr ] = trans_sys.rollout( y0, z0, g, scl * input_arr, t0i, t_arr  );
    
    % The original demonstrated trajectory
    if i == 1 
        plot3( a, y_arr( 1, : ), y_arr( 2, : ), y_arr( 3, : ), 'linewidth', 5, 'color', [0.8500 0.3250 0.0980] )
    else
        plot3( a, y_arr( 1, : ), y_arr( 2, : ), y_arr( 3, : ), 'linewidth', 5, 'color', 'black' )
    end

    scatter3( a, g( 1 ), g( 2 ), g( 3 ), 200, 'filled', 'markeredgecolor', [0 0.4470 0.7410], 'markerfacecolor', 'white', 'linewidth', 4 );

end


%%  -- (1F) Saving the Data, the weights for the learned trajectory

data = struct;

% Parameters of DMP and Learned Weighted
data.name    = 'lift_up_down';
data.tau     = tau;
data.alpha_s = alpha_s;
data.alpha_z = alpha_z;
data.beta_z  =  beta_z;
data.weight  = w_arr;

% Initial and Goal Location
data.goal =  p_arr_filt( :, end );
data.z0   = dp_arr_filt( :, 1   )/tau;

save( [ 'learned_parameters/', traj_names{ idx } , '_pos.mat' ], 'data' );