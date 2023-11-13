% [Title]     Data Filter and learn the Weight Array
% [Author]    Moses Chong-ook Nah
% [Email]     mosesnah@mit.edu
% [Update]    At 2023.11.12

%% (--) Initialization
clear; close all; clc;
fig_config( 'fontSize', 20, 'markerSize', 10 )

%% =======================================================
%% (1-) Position Data Filtering and Learning the Weights
%%  -- (1A) Calling the data and get the Forward Kinematics 

file_name = './data/example4/iiwa_example_pos';
raw_data = parse_txt( [ file_name, '.txt' ], 0 );
t_arr_demo = raw_data( :, 1 )' - raw_data( 1, 1 );
q_arr_demo = raw_data( :, 2:8 )';

% Call the Forward Kinematics of the Robot
robot = iiwa14( 'high' );
robot.init( );

[ ~, Np ] = size( q_arr_demo );

% End-effect Forward Kinematics, both position and orientation
% We will learn the position, but also saving the orientation for backup
p_arr_raw = zeros( 3, Np );
R_arr_raw = zeros( 3, 3, Np );

for i = 1 : Np
    H = robot.getForwardKinematics( q_arr_demo( :, i ) );
    p_arr_raw( :, i )    = H( 1:3,   4 );
    R_arr_raw( :, :, i ) = H( 1:3, 1:3 ); 
end

is_check = true;

if is_check
    f = figure( ); a = axes( 'parent', f );
    hold( a, 'on' ); axis equal
    plot3( a, p_arr_raw( 1, : ), p_arr_raw( 2, : ), p_arr_raw( 3, : ), 'linewidth', 3 )
end

%%  -- (1B) Filtering out the Position data

% Simple time diff
dp_arr_raw  = data_diff(  p_arr_raw );
ddp_arr_raw = data_diff( dp_arr_raw );

% Filtering out the Data
  p_arr_filt = zeros( size( p_arr_raw ), 'like', p_arr_raw );
 dp_arr_filt = zeros( size( p_arr_raw ), 'like', p_arr_raw );
ddp_arr_filt = zeros( size( p_arr_raw ), 'like', p_arr_raw );

% Filtered position, velocity and acceleration
for i = 1 : 3
     p_arr_filt( i, : ) = smoothdata( p_arr_raw( i, : ), "gaussian", 50 );

   tmp_vel = data_diff( p_arr_filt );
    dp_arr_filt( i, : ) = smoothdata( tmp_vel( i, : ), "gaussian", 50 );

   tmp_acc = data_diff( dp_arr_filt );
   ddp_arr_filt( i, : ) = smoothdata( tmp_acc( i, : ), "gaussian", 50 );   
end

is_check = false;

if is_check

    tmp_plot = cell( 1, 6 );
    tmp_plot{ 1 } = p_arr_raw;
    tmp_plot{ 2 } = p_arr_filt;
    
    tmp_plot{ 3 } = dp_arr_raw;
    tmp_plot{ 4 } = dp_arr_filt;
    
    tmp_plot{ 5 } = ddp_arr_raw;
    tmp_plot{ 6 } = ddp_arr_filt;
    
    tmp_label = [ "Pos.", "Vel.", "Acc." ];
    
    % Row
    for i = 1:3

        % Column
        for j = 1 :2
            subplot( 3, 2, 2*(i-1)+j )
            hold on
            plot( t_arr, cell{ 2*(i-1)+j }, 'linewidth', 3 )
            set( gca, 'xlim', [ 0, max( t_arr ) ], 'xticklabel', {}, 'fontsize', 30 )
            
            if j == 1
               ylabel( tmp_label( i ) )
               tmp = get( gca, 'ylim' );
            else
               set( gca, 'ylim', tmp ) 
            end
        end
    
    end

end    

%%  -- (1C) Learning the Weights and Double Check

% Parameters of the DMP
alpha_s   =  2.0;
alpha_z   = 2000.0;
beta_z    = 0.5 * alpha_z;
N         = 100;
tau       = max( t_arr_demo );

% The Three elements of DMP
cs        = CanonicalSystem( 'discrete', tau, alpha_s );
trans_sys = TransformationSystem( alpha_z, beta_z, cs );
fs        = NonlinearForcingTerm( cs, N );

% Calculating the required Nonlinear Forcing Term
P  = length( t_arr_demo );
y0d = p_arr_filt( :, 1   );
gd  = p_arr_filt( :, end );
f_arr = trans_sys.get_desired( p_arr_filt, dp_arr_filt, ddp_arr_filt, gd );

% The phi matrix 
Phi_mat = zeros( P, N );

% Interating along the sample points
for i = 1 : P 
    t = t_arr_demo( i );
    Phi_mat( i, : ) = fs.calc_multiple_ith( t, 1:N )/ fs.calc_whole_at_t( t ) * cs.calc( t );
end

% Linear Least-square fitting
w_arr = transpose( ( Phi_mat' * Phi_mat )^(-1) * Phi_mat' * f_arr' );

% Double check with the learned weights
is_check = true;

if is_check 
    y0 =  p_arr_filt( :,   1 );
    z0 = dp_arr_filt( :,   1 )/tau;
    g  =  p_arr_filt( :, end );
    t0i = 0.0;
    T   = max( t_arr_demo );
    
    input_arr = fs.calc_forcing_term( t_arr_demo( 1:end-1 ), w_arr, t0i, eye( length( g ) ) );
    
    [ y_arr, z_arr, dy_arr ] = trans_sys.rollout( y0, z0, g, input_arr, t0i, t_arr_demo  );
    
    f = figure( ); a = axes( 'parent', f );
    hold( a, 'on' );
    plot( a, t_arr_demo, p_arr_filt, 'linewidth', 10, 'linestyle', '--', 'color', 'black' )
    plot( a, t_arr_demo,      y_arr, 'linewidth',  5, 'linestyle',  '-', 'color', 'blue'  )
end

%%  -- (1D) Generalize the demonstrated trajectory 

% Since the weights are learned, we can scale and rotate the trajectory 
% Rotation
y0 = zeros( 3, 1 );
g  =  p_arr_filt( :, end ) - p_arr_filt( :, 1 ); 
z0 = dp_arr_filt( :, 1 )/cs.tau;

t0i   = 0.0;
T     = 5.0;
t_arr = 0:1e-3:T;

f = figure( ); a = axes( 'parent', f );
hold( a, 'on' ); axis equal
set( a, 'view', [46.1836, 38.4642], 'visible', 'off' )
scatter3( a, 0, 0, 0, 200, 'filled', 'markeredgecolor', 'black',  'markerfacecolor', 'white', 'linewidth', 2 );

theta_arr = 1:60:360;
scl_arr   = [ 1.0, 2.0, 4.0 ];
% for i = 1 : length( theta_arr )
for i = 1:length( scl_arr )

    % theta = theta_arr( i );
    scl = scl_arr( i );

    Rmat = rotz( 0 );
        
    g  = p_arr_filt( :, end ) - p_arr_filt( :, 1 ); 
    g  = scl * Rmat * g;
     
    input_arr = fs.calc_forcing_term( t_arr( 1:end-1), w_arr, t0i, eye( length( g ) ) );
    
    [ y_arr, z_arr, dy_arr ] = trans_sys.rollout( y0, z0, g, scl * Rmat * input_arr, t0i, t_arr  );
    
    if i == 1
        plot3( a, y_arr( 1, : ), y_arr( 2, : ), y_arr( 3, : ), 'linewidth', 5, 'color', [0.8500 0.3250 0.0980] )
    else
        plot3( a, y_arr( 1, : ), y_arr( 2, : ), y_arr( 3, : ), 'linewidth', 5, 'color', 'black' )
    end

    scatter3( a, g( 1 ), g( 2 ), g( 3 ), 200, 'filled', 'markeredgecolor', [0 0.4470 0.7410], 'markerfacecolor', 'white', 'linewidth', 4 );

end

%% =======================================================
%% (2-) SO3 Data Filtering and Learning the Weights
%%  -- (2A) Calling the data and get the Forward Kinematics 

clear all; close all; clc;
file_name = './data/example3/iiwa_example_orient';
raw_data = parse_txt( [ file_name, '.txt' ], 0 );
t_arr_raw = raw_data( :, 1 )' - raw_data( 1, 1 );
q_arr_raw = raw_data( :, 2:8 )';

% Call the robot data for the Forward Kinematics Map
robot = iiwa14( 'high' );
robot.init( );

[ ~, Np ] = size( q_arr_raw );

% End-effect Forward Kinematics, position
p_arr_raw  = zeros( 3, Np );
p_arr_filt = zeros( 3, Np );

% End-effect Forward Kinematics, orientation
R_arr_raw  = zeros( 3, 3, Np );
R_arr_filt = zeros( 3, 3, Np );

% Get the filtered q and dq values
q_arr_filt  = zeros( size( q_arr_raw ), 'like', q_arr_raw );
dq_arr_filt = zeros( size( q_arr_raw ), 'like', q_arr_raw );

% Over the degrees of freedom
for i = 1 : 7
     q_arr_filt( i, : ) = smoothdata(  q_arr_raw( i, : ), "gaussian", 50 );

    dq_arr_filt( i, : ) = data_diff( q_arr_filt( i , : ) ); 
    dq_arr_filt( i, : ) = smoothdata( dq_arr_filt( i, : ), "gaussian", 50 );
end

% Get the Forward Kinematics Map
for i = 1 : Np
    H = robot.getForwardKinematics( q_arr_raw( :, i ) );
    R_arr_raw( :, :, i ) = H( 1:3, 1:3 ); 
    p_arr_raw( :, i ) = H( 1:3,   4 ); 

    H = robot.getForwardKinematics( q_arr_filt( :, i ) );
    R_arr_filt( :, :, i ) = H( 1:3, 1:3 );
    p_arr_filt( :, i ) = H( 1:3,   4 ); 
end

% Getting the body angular velocity 
w0_arr_filt = zeros( 3, Np );

for i = 1 : Np
    JB  = robot.getBodyJacobian( q_arr_filt( :, i ), 7 );
    JBr = JB( 4:6, : );
    w0_arr_filt( :, i ) = JBr * dq_arr_filt( :, i );
end

% Getting the body acceleration
dw0_arr_filt = zeros( 3, Np );

for i = 1 : 3
    dw0_arr_filt(  i, : ) =  data_diff( w0_arr_filt( i, : ) ); 
    dw0_arr_filt( i, : ) = smoothdata( dw0_arr_filt( i, : ), "gaussian", 50 );
end

is_check = false;

if is_check
    f = figure( ); a = axes( 'parent', f );
    hold( a, 'on' ); axis equal
    set( a, 'view',  [69.2886, 5.1376] )
    plot3( a, p_arr_raw( 1, : ), p_arr_raw( 2, : ), p_arr_raw( 3, : ), 'linewidth', 3, 'color', 'black' )

    scl = 0.02;
    for i = 1 : 100 : Np

        p = p_arr_raw( :, i );
        R = scl * R_arr_raw( :, :, i );
        
        scatter3( a, p( 1 ), p( 2 ), p( 3 ), 500, 'filled', 'markerfacecolor', 'white', 'markeredgecolor', 'black', 'linewidth', 5)

        quiver3( a, p( 1 ), p( 2 ), p( 3 ), R( 1, 1 ), R( 2, 1 ), R( 3, 1 ), 'linewidth',3, 'color', 'r')
        quiver3( a, p( 1 ), p( 2 ), p( 3 ), R( 1, 2 ), R( 2, 2 ), R( 3, 2 ), 'linewidth',3, 'color', 'g')
        quiver3( a, p( 1 ), p( 2 ), p( 3 ), R( 1, 3 ), R( 2, 3 ), R( 3, 3 ), 'linewidth',3, 'color', 'b')

        p = p_arr_filt( :, i );
        R = scl * R_arr_filt( :, :, i );
        
        scatter3( a, p( 1 ), p( 2 ), p( 3 ), 500, 'filled', 'markerfacecolor', 'white', 'markeredgecolor', 'black', 'linewidth', 5)

        quiver3( a, p( 1 ), p( 2 ), p( 3 ), R( 1, 1 ), R( 2, 1 ), R( 3, 1 ), 'linewidth',3, 'color', 'r')
        quiver3( a, p( 1 ), p( 2 ), p( 3 ), R( 1, 2 ), R( 2, 2 ), R( 3, 2 ), 'linewidth',3, 'color', 'g')
        quiver3( a, p( 1 ), p( 2 ), p( 3 ), R( 1, 3 ), R( 2, 3 ), R( 3, 3 ), 'linewidth',3, 'color', 'b')

    end

end

%%  -- (2B) Learning the Weights

% Get the goal rotation matrix 
Rg = R_arr_filt( :, :, end );
R_arr_del = zeros( size( R_arr_filt ), 'like', R_arr_filt );

% Getting the delta of SO3
for i = 1 : Np
    R_arr_del( :, :, i ) = R_arr_filt( :, :, i )' * Rg;
end

% Define the transformation system and others to learn the weights

% Parameters of the DMP
alpha_s   = 2.0;
alpha_z   = 2000.0;
beta_z    = 0.5 * alpha_z;
N         = 100;
tau       = max( t_arr_raw );

% The Three elements of DMP
cs        = CanonicalSystem( 'discrete', tau, alpha_s );
trans_sys = TransformationSystemSO3( alpha_z, beta_z, cs );
fs        = NonlinearForcingTerm( cs, N );

f_arr = trans_sys.get_desired( R_arr_filt, w0_arr_filt, dw0_arr_filt );


%% =======================================================
%% (3-) Angular Velocity Data Filtering and Learning the Weights
%% Filter out the w_arr, dw_arr 

for i = 1 : N
    tmp = robot.getHybridJacobian( q_arr( :, i ) );
    JH_arr( :, :, i )  = tmp( 4:6, : );
end

% The raw w_arr 
 w_arr_raw  = zeros( size( p_arr_raw ), 'like', p_arr_raw );
dw_arr_raw  = zeros( size( p_arr_raw ), 'like', p_arr_raw );
 w_arr_filt = zeros( size( p_arr_raw ), 'like', p_arr_raw );
dw_arr_filt = zeros( size( p_arr_raw ), 'like', p_arr_raw );

% The non-filtered dq_arr
dq_arr_raw  = data_diff( q_arr );
dq_arr_filt = zeros( size( q_arr ), 'like', q_arr );

for i = 1 : 7
    dq_arr_filt( i, : ) = smoothdata( dq_arr_raw( i, : ), "gaussian", 50 );
end

for i = 1 : N
     w_arr_raw( :, i ) = JH_arr( :, :, i ) * dq_arr_raw( :, i );
    w_arr_filt( :, i ) = JH_arr( :, :, i ) * dq_arr_filt( :, i );
end

 dw_arr_raw = data_diff( w_arr_raw  );
dw_arr_filt = data_diff( w_arr_filt );

for i = 1 : 3
    dw_arr_filt( i, : ) = smoothdata( dw_arr_filt( i, : ), "gaussian", 50 );
end
f = figure( );
subplot( 2, 2, 1 )
plot( t_arr, w_arr_raw, 'linewidth', 3 )
set( gca, 'xlim', [0, max( t_arr ) ] )
ylim1 = get( gca, 'ylim' );
ylabel( 'Ang. Vel. ($rad/s$)')

subplot( 2, 2, 2 )
plot( t_arr, w_arr_filt, 'linewidth', 3 )
set( gca, 'xlim', [0, max( t_arr ) ], 'ylim', ylim1 )

subplot( 2, 2, 3 )
plot( t_arr, dw_arr_raw, 'linewidth', 3 )
set( gca, 'xlim', [0, max( t_arr ) ], 'ytick', [-5, 0, 5, 10]*10^-5 )
ylim2 = get( gca, 'ylim' );
ylabel( 'Ang. Acc. ($rad/s^2$)')
xlabel( 'Time ($s$)')

subplot( 2, 2, 4 )
plot( t_arr, dw_arr_filt, 'linewidth', 3 )
set( gca, 'xlim', [0, max( t_arr ) ], 'ylim', ylim2, 'ytick', [-5, 0, 5, 10]*10^-5 )
xlabel( 'Time ($s$)')

%% Check for sign flip of the quaternion

tmp = diff( quat_arr_raw' )';

% If there exists a value that changes so drastically, find the indexes 
[r, c_arr, ] = find( abs( tmp ) >= 1 );

% If not empty
if( ~isempty( c_arr ) ) 
    % c is the index, flip the sign 
    for c = c_arr
        quat_arr_raw( :,c+1:end ) = -quat_arr_raw( :,c+1:end );
    end
end

%% Save/Export Data

% Saving the crucial data 
save( [file_name, '.mat'], 'w_arr_raw', 'w_arr_filt', 'dw_arr_raw',  'dw_arr_filt', ...
                           'p_arr_raw', 'p_arr_filt', 'dp_arr_raw', 'dp_arr_filt', 'ddp_arr_raw', 'ddp_arr_filt', ... 
                           'R_arr_raw', 'quat_arr_raw', 'q_arr', 't_arr' )
