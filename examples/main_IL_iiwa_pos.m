%% ==================================================================
%% [Title] Imitation Learning, in 3D Space, Position
% Author: Moses Chong-ook Nah
%  Email: mosesnah@mit.edu
%   Date: 2023.08.15
%% ==================================================================

%% [0A] Initialization
clear; close all; clc;
fig_config( 'fontSize', 20, 'markerSize', 10 )

%% ---- [1A] Call iiwa14 using Explicit-MATLAB

% Set figure size and attach robot to simulation
robot = iiwa14( 'high' );
robot.init( );

%% ---- [1B] Load the Data we aim to imitate

dir_name  = './data/example6/';
file_name = 'test3';
data_raw = load( [ dir_name, file_name, '.mat' ] );


%% ---- [1C] Imitation Learning - Learning the Weights

% Shorten the name of the variables
t_arr = data_raw.t_arr;

ddp_des = data_raw.ddp_arr_filt;
dp_des  =  data_raw.dp_arr_filt;
p_des   =   data_raw.p_arr_filt;

N       = 100;              % Number of Basis Function
alpha_z = 1000.0;            % Gain of alpha_z
beta_z  = 1/4 * alpha_z;    % Gain of  beta_z
alpha_s = 1.0;              % Gain of Canonical System
tau     = max( t_arr );     % Time array
y0      =  p_des( :, 1 );   % Initial Position
z0      = dp_des( :, 1 );   % Initial Velocity

% Get the number of total samples
[ ~, P_total ] = size( p_des );

% Don't need to use all the data
nP_step = 3;
idx = 1 : nP_step : P_total;

% The number of data points for Imitation Learning
P = length( idx );

% Get the goal location as the final position
g = p_des( :, idx( end ) );

% Create a single Canonical System and Three Nonlinear Forcing Terms
cs     = CanonicalSystem( 'discrete', tau, alpha_s );
fs_x   = NonlinearForcingTerm( cs, N );
fs_y   = NonlinearForcingTerm( cs, N );
fs_z   = NonlinearForcingTerm( cs, N );

tsys_x = TransformationSystem( alpha_z, beta_z, tau, y0( 1 ), z0( 1 ) );
tsys_y = TransformationSystem( alpha_z, beta_z, tau, y0( 2 ), z0( 2 ) );
tsys_z = TransformationSystem( alpha_z, beta_z, tau, y0( 3 ), z0( 3 ) );


% Save the Nonlinear Forcing Terms as a cell
fs = cell( 1, 3 ); tsys = cell( 1, 3 );
  fs{ 1 } =   fs_x;   fs{ 2 } =   fs_y;   fs{ 3 } =   fs_z;
tsys{ 1 } = tsys_x; tsys{ 2 } = tsys_y; tsys{ 3 } = tsys_z;

% Learning the N weights based on Locally Weighted Regression
w_arr_LWR = zeros( 3, N );

% Iterating over the three dimension and array
% For each Degree of freedom
t_arr_P = t_arr( idx ); 
for k = 1 : 3
    a_arr = ( g( k ) - y0( k ) ) * cs.calc( t_arr_P );
    f_arr = tsys{ k }.get_desired( p_des( k, idx ), dp_des( k, idx ), ddp_des( k, idx ), g( k ) );
    
    for i = 1: N
        phi_arr = fs{ k }.calc_ith( t_arr_P, i );    
        w_arr_LWR( k, i ) = sum( a_arr .* phi_arr .* f_arr ) / sum( a_arr .* phi_arr .* a_arr );
    end
end

% Also Learning the weights from Least-square Solution
w_arr_LSS = zeros( 3, N );

for k = 1 : 3

    f_arr = tsys{ k }.get_desired( p_des( k, idx ), dp_des( k, idx ), ddp_des( k, idx ), g( k ) );
    phi_mat = zeros( P, N );

    for i = 1 : P
        phi_sum = fs{ k }.calc_whole_at_t( t_arr_P( i ) );
        for j = 1 : N
            phi_mat( i, j ) = fs{ k }.calc_ith( t_arr_P( i ), j ) / phi_sum * ( g( k ) - y0( k ) ) * cs.calc( t_arr_P( i ) );
        end
    end

    % Get w_arr with Least square solution
    tmp = ( phi_mat' * phi_mat )^(-1) * phi_mat' * f_arr';
    w_arr_LSS( k, : ) = tmp';    
end

%% ---- [1D] Generating a Full trajectory with the Transformation System

% The time step of the simulation and its number of iteration
dt = t_arr( 2) - t_arr( 1 );
Nt = P_total * 2;

% For plotting and saving the data
t_arr2 = dt *( 0:(Nt-1) );
y_arr  = zeros( 3, Nt );
z_arr  = zeros( 3, Nt );
dy_arr = zeros( 3, Nt );
dz_arr = zeros( 3, Nt );

t = 0;
t0i = 3.0;

tsys{ 1 }.reset( );
tsys{ 2 }.reset( );
tsys{ 3 }.reset( );

scl = 0.8;

for i = 0 : (Nt-1)
    
    if t <= t0i
        for k = 1 : 3
            y_arr( k, i + 1 )  =  y0( k );
            z_arr( k, i + 1 )  =  z0( k );
        end

    elseif t0i <= t && t <= t0i + max( t_arr )

        ttmp = t - t0i;

        for k = 1 : 3
            % Calculating the input from the weights
            % First, check if whole activation value is 0
            phi_sum = fs{ k }.calc_whole_at_t( ttmp );
    
            if phi_sum ~= 0
                f_input = fs{ k }.calc_whole_weighted_at_t( ttmp, w_arr_LSS( k, : ) )/phi_sum;
                f_input = f_input*( g( k )-y0( k ) )*cs.calc( ttmp );
            end

     
            [ y, z, dy, dz ] = tsys{ k }.step( g( k ), f_input, dt );
            y_arr(  k, i + 1 )  =  y;
            z_arr(  k, i + 1 )  =  z;
            dy_arr(  k, i + 1 ) = dy;
            dz_arr( k, i + 1 )  = dz;    
        end

    else
        for k = 1 : 3
            [ y, z, dy, dz ] = tsys{ k }.step( g( k ), 0, dt );
            y_arr(  k, i + 1 )  =  y;
            z_arr(  k, i + 1 )  =  z;
            dy_arr(  k, i + 1 ) = dy;
            dz_arr( k, i + 1 )  = dz;    
        end
    end
    t = t + dt;
end

a1 = subplot( 3, 2, 1 );
hold( a1, 'on' );
plot( a1, t_arr2, y_arr( 1, : ), 'linewidth', 5, 'color', [0.0000, 0.4470, 0.7410] )
plot( a1, t_arr, data_raw.p_arr_filt( 1, : ), 'linewidth', 3, 'linestyle', '--', 'color', 'k' )
set( a1, 'xlim', [0, max( t_arr2 )], 'fontsize', 25, 'xticklabel', {} )
title( a1, 'X-Position', 'Fontsize', 30)

a2 = subplot( 3, 2, 3 );
hold( a2, 'on' );
plot( a2, t_arr2, y_arr( 2, : ), 'linewidth', 5, 'color', [0.8500, 0.3250, 0.0980] )
plot( a2, t_arr, data_raw.p_arr_filt( 2, : ), 'linewidth', 3, 'linestyle', '--', 'color', 'k' )
set( a2, 'xlim', [0, max( t_arr2 )],  'fontsize', 25, 'xticklabel', {} )
title( a2, 'Y-Position', 'Fontsize', 30);

a3 = subplot( 3, 2, 5 );
hold( a3, 'on' );
plot( a3, t_arr2, y_arr( 3, : ), 'linewidth', 5, 'color', [0.9290 0.6940 0.1250] )
plot( a3, t_arr, data_raw.p_arr_filt( 3, : ), 'linewidth', 3, 'linestyle', '--', 'color', 'k' )
set( a3, 'xlim', [0, max( t_arr2 )], 'fontsize', 25 )
xlabel( a3, 'Time (sec)', 'Fontsize', 30)
title( a3, 'Z-Position', 'Fontsize', 30);


dy_arr = data_diff( y_arr );

a4 = subplot( 3, 2, 2 );
hold( a4, 'on' );
plot( a4, t_arr2, dy_arr( 1, : ), 'linewidth', 5, 'color', [0.0000, 0.4470, 0.7410] )
plot( a4, t_arr, data_raw.dp_arr_filt( 1, : ), 'linewidth', 3, 'linestyle', '--', 'color', 'k' )
set( a4, 'xlim', [0, max( t_arr2 )], 'fontsize', 25, 'xticklabel', {} )
title( a4, 'X-Velocity', 'Fontsize', 30)

a5 = subplot( 3, 2, 4 );
hold( a5, 'on' );
plot( a5, t_arr2, dy_arr( 2, : ), 'linewidth', 5, 'color', [0.8500, 0.3250, 0.0980] )
plot( a5, t_arr, data_raw.dp_arr_filt( 2, : ), 'linewidth', 3, 'linestyle', '--', 'color', 'k' )
set( a5, 'xlim', [0, max( t_arr2 )],  'fontsize', 25, 'xticklabel', {} )
title( a5, 'Y-Velocity', 'Fontsize', 30);

a6 = subplot( 3, 2, 6 );
hold( a6, 'on' );
plot( a6, t_arr2, dy_arr( 3, : ), 'linewidth', 5, 'color', [0.9290 0.6940 0.1250] )
plot( a6, t_arr, data_raw.dp_arr_filt( 3, : ), 'linewidth', 3, 'linestyle', '--', 'color', 'k' )
set( a6, 'xlim', [0, max( t_arr2 )], 'fontsize', 25 )
xlabel( a6, 'Time (sec)', 'Fontsize', 30)
title( a6, 'Z-Velocity', 'Fontsize', 30);



%% Compare in Explicit

anim = Animation( 'Dimension', 3, 'xLim', [-0.7,0.7], 'yLim', [-0.7,0.7], 'zLim', [0,1.4], 'isSaveVideo', true );
anim.init( );
anim.attachRobot( robot )  
    
plot3(    anim.hAxes, p_des( 1, :   ), p_des( 2, :   ), p_des( 3, :   ), 'linewidth', 3, 'linestyle', '--', 'color', 'k' );
% plot3(    anim.hAxes, y_arr( 1, :   ), y_arr( 2, :   ), y_arr( 3, :   ), 'linewidth', 3, 'linestyle', '-', 'color', [0.4940 0.1840 0.5560]    );
scatter3( anim.hAxes, p_des( 1, 1   ), p_des( 2, 1   ), p_des( 3, 1   ), 300, 'filled', 'linewidth', 3, 'markerfacecolor', [0 0.4470 0.7410], 'markeredgecolor', 'k'  );
scatter3( anim.hAxes, p_des( 1, end ), p_des( 2, end ), p_des( 3, end ), 300, 'filled', 'linewidth', 3, 'markerfacecolor', [0 0.4470 0.7410], 'markeredgecolor', 'k' );

view([90,0])

q_arr = data_raw.q_arr; 

for i = 1 : length( t_arr )
    robot.updateKinematics( q_arr( :, i ) );
    H = robot.getForwardKinematics( q_arr( :, i ) );
    anim.update( t_arr( i ) );
end

anim.close( );


%% If ready, generate data input for robot

% At the end-of the day, what is matters is the delta movement
% Saving both the position and velocity values 
writematrix(  y_arr - y_arr( :, 1 ), [ dir_name, 'pos_data.csv' ] ) 

% Saving the velocity data too
writematrix( dy_arr, [ dir_name, 'vel_data.csv' ] ) 
