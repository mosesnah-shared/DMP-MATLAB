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

dir_name  = './data/example3/';
file_name = 'iiwa_example_orient';
data_raw = load( [ dir_name, file_name, '.mat' ] );

% Array of quaternion, angular velocity and acceleration
t_arr      = data_raw.t_arr; 
quat       = data_raw.quat_arr_raw;
omega_arr  = data_raw.w_arr_filt;
domega_arr = data_raw.dw_arr_filt;


% For Imitation Learning, one should define the number of basis function
N  = 30;

% Parameters of the 3 DMPs
alpha_z = 400.0;
alpha_s = 1.0;
beta_z  = 600.0;
tau     = max( t_arr );

% Get the number of total samples
[ ~, P_total ] = size( quat );

idx = 1 : 5 : P_total;
P = length( idx );

% Get the goal location  
g     = quat( :, idx( end ) );
cs    = CanonicalSystem( 'discrete', tau, alpha_s );
fs_w  = NonlinearForcingTerm( cs, N );
% Learning the N weights based on Locally Weighted Regression
w_arr = zeros( 3, N );

vec0 = quat_log( quat_mul( g, quat_conj( quat( :, 1 ) ) ) ); 

phi_arr = zeros( 3*P, 3*N );
f_arr   = zeros( 3*P, 1 );

for j = 1 : P

    ddy_des = domega_arr( :, idx( j ) );
     dy_des =  omega_arr( :, idx( j ) ); 
    quat_err = quat_log( quat_mul( g, quat_conj( quat( :, idx( j ) ) ) ) ); 

    f_arr( 3*(j-1)+1:3*j ) = tau * ddy_des + alpha_z * dy_des - alpha_z * beta_z * quat_err;

end

for i = 1 : N
    for j = 1 : P
        
        % Get the sum of phis
        phi_act_arr = zeros( 1, N );
        for k = 1 : N
            phi_act_arr( k ) = fs_w.calc_ith(  t_arr( idx( j ) ), k );
        end

        phi_arr( 3*(j-1)+1:3*j, 3*(i-1)+1:3*i) = fs_w.calc_ith(  t_arr( idx( j ) ), i ) / sum( phi_act_arr ) * cs.calc( t_arr( idx( j ) ) ) * eye( 3 );
    end
end

w_arr = inv( phi_arr' * phi_arr ) * phi_arr' * f_arr;

%% Using the weights for the transformation system

trans = TransformationSystem_quat(  alpha_z, beta_z, tau, quat( :, 1 ), omega_arr( :, 1 ) );

% The time step of the simulation and its number of iteration
dt = t_arr( 2 ) - t_arr( 1 );
Nt = P_total * 3;

% The total time and its time array
T  = dt * Nt;
t_arr = dt * (0:(Nt-1));

% For plotting the orientation 
quat_data_arr  = zeros( 4, Nt );
omega_data_arr = zeros( 3, Nt );

t = 0;

f_input_arr = zeros( 1, Nt );

y0 = quat( :, 1 );
z0 = omega_arr( :, 1 );

y = y0;
z = z0;

t0i = 3.0;

for i = 1 : Nt
    
    if t <= t0i

        quat_data_arr(  :, i  ) = y0;
        omega_data_arr( :, i  ) = z0;

    else    
        % Calculating the input from the weights        
        vectmp = quat_log( quat_mul( g, quat_conj( y ) ) ); 
        
        ttmp = t - t0i;

        % Get f_input 
        % Get the sum of phis
        phi_act_arr = zeros( 1, N );    
        for k = 1 : N
            phi_act_arr( k ) = fs_w.calc_ith( ttmp, k );
        end
        
        f_input = zeros( 3, 1 );
        for k = 1 : N
            f_input = f_input + fs_w.calc_ith( ttmp, k ) * w_arr( 3*(k-1)+1:3*k ) / sum( phi_act_arr ) * cs.calc( ttmp );
        end
    
        [ quat_new, w_new, ~, ~ ] = trans.step( g, f_input, dt );
        quat_data_arr(  :, i  ) = quat_new;
        omega_data_arr( :, i  ) = w_new;

        y = quat_new;
    end

    
    t = t + dt;
end


%% Change quat to SO3

[ ~, Nr] = size( quat_data_arr );

R_data_arr = zeros( 3, 3, Nr );

for i = 1 : Nr
    R_data_arr( :, :, i ) = quat_to_SO3( quat_data_arr( :, i ) );
end

c_arr = [0.0000, 0.4470, 0.7410;
         0.8500, 0.3250, 0.0980;
         0.9290, 0.6940, 0.1250;
         0.4940  0.1840, 0.5560];

ylabels = { '$q_w$', '$q_x$', '$q_y$', '$q_z$'};

for i = 1 : 4
    subplot( 4, 2, (2*i-1) )
    hold on
    plot( t_arr , quat_data_arr( i, :), 'color', 'k', 'linewidth', 3 );
    plot( data_raw.t_arr, quat( i, :) , 'linestyle', '--', 'color', c_arr( i, : ), 'linewidth', 5);
    set( gca, 'xlim', [ 0, max( t_arr ) ], 'fontsize', 25, 'ylim', [-1,1] )
    if i ~= 4
        set( gca, 'xticklabel', {})
    end

    if i == 4
        xlabel( 'Time (sec)', 'fontsize', 30 )
    end
    ylabel( ylabels{ i }, 'fontsize', 30 )
end

%% Compare in Explicit

anim = Animation( 'Dimension', 3, 'xLim', [-0.7,0.7], 'yLim', [-0.7,0.7], 'zLim', [0,1.4], 'isSaveVideo', false );
anim.init( );
anim.attachRobot( robot )  

q_arr = data_raw.q_arr;
[ ~, N ] = size( q_arr );

p_arr     = zeros( 3, N );
R_arr     = zeros( 3, 3, N );
R_arr_del = zeros( 3, 3, N );

for i = 1 : N
    robot.updateKinematics( q_arr( :, i ) );
    H = robot.getForwardKinematics( q_arr( :, i ) );
    p_arr( :, i ) = H( 1:3, 4 );
    R_arr( :, :, i ) = H( 1:3, 1:3 );

end

ttmp1 = R_arr( :, :, 200 );
for i = 1 : N
    R_arr_del( :, :, i ) = R_arr( :, :, 1 )' * R_arr( :, :, i );
end

view([90,0])

q_arr = data_raw.q_arr; 

plot3( anim.hAxes,  p_arr( 1, : ), p_arr( 2, : ), p_arr( 3, : ), 'linewidth', 3, 'linestyle', '--', 'color', 'k')
scl = 0.1;

for i = 1:200:N
    quiver3( anim.hAxes,  p_arr( 1, i ), p_arr( 2, i ), p_arr( 3, i ), scl * R_arr( 1, 1, i ), scl * R_arr( 2, 1, i ), scl * R_arr( 3, 1, i ), 'linewidth', 3, 'color', 'r')
    quiver3( anim.hAxes,  p_arr( 1, i ), p_arr( 2, i ), p_arr( 3, i ), scl * R_arr( 1, 2, i ), scl * R_arr( 2, 2, i ), scl * R_arr( 3, 2, i ), 'linewidth', 3, 'color', 'g')
    quiver3( anim.hAxes,  p_arr( 1, i ), p_arr( 2, i ), p_arr( 3, i ), scl * R_arr( 1, 3, i ), scl * R_arr( 2, 3, i ), scl * R_arr( 3, 3, i ), 'linewidth', 3, 'color', 'b')
end

offset = 0.25;
plot3( anim.hAxes,  p_arr( 1, : ), p_arr( 2, : ), offset +  p_arr( 3, : ), 'linewidth', 3, 'linestyle', '--', 'color', 'k')


for i = 1:200:N
    R_tmp = quat2rotm( quat_data_arr( :, i )' );
    quiver3( anim.hAxes,  p_arr( 1, i ), p_arr( 2, i ), offset + p_arr( 3, i ), scl * R_arr_del( 1, 1, i ), scl * R_arr_del( 2, 1, i ), scl * R_arr_del( 3, 1, i ), 'linewidth', 3, 'color', 'r')
    quiver3( anim.hAxes,  p_arr( 1, i ), p_arr( 2, i ), offset + p_arr( 3, i ), scl * R_arr_del( 1, 2, i ), scl * R_arr_del( 2, 2, i ), scl * R_arr_del( 3, 2, i ), 'linewidth', 3, 'color', 'g')
    quiver3( anim.hAxes,  p_arr( 1, i ), p_arr( 2, i ), offset + p_arr( 3, i ), scl * R_arr_del( 1, 3, i ), scl * R_arr_del( 2, 3, i ), scl * R_arr_del( 3, 3, i ), 'linewidth', 3, 'color', 'b')
end

offset = 0.5;
plot3( anim.hAxes,  p_arr( 1, : ), p_arr( 2, : ),offset +  p_arr( 3, : ), 'linewidth', 3, 'linestyle', '--', 'color', 'k')

for i = 1:200:N
    R_tmp = quat2rotm( quat_data_arr( :, i )' );
    quiver3( anim.hAxes,  p_arr( 1, i ), p_arr( 2, i ), offset + p_arr( 3, i ), scl * R_tmp( 1, 1 ), scl * R_tmp( 2, 1 ), scl * R_tmp( 3, 1 ), 'linewidth', 3, 'color', 'r')
    quiver3( anim.hAxes,  p_arr( 1, i ), p_arr( 2, i ), offset + p_arr( 3, i ), scl * R_tmp( 1, 2 ), scl * R_tmp( 2, 2 ), scl * R_tmp( 3, 2 ), 'linewidth', 3, 'color', 'g')
    quiver3( anim.hAxes,  p_arr( 1, i ), p_arr( 2, i ), offset + p_arr( 3, i ), scl * R_tmp( 1, 3 ), scl * R_tmp( 2, 3 ), scl * R_tmp( 3, 3 ), 'linewidth', 3, 'color', 'b')
end

for i = 1 : N
    robot.updateKinematics( q_arr( :, i ) );
    H = robot.getForwardKinematics( q_arr( :, i ) );
    anim.update( t_arr( i ) );
end

anim.close( );


%% Final test of R_arr_del 



%% If ready, generate Save data output for robot.

% Reshaping the R_arr as 2D array
[ ~, ~, N ] = size( R_data_arr );
R_arr_del = zeros( 3, 3, N );

for i = 1 : N
    R_arr_del( :, :, i ) = R_data_arr( :, :, 1 )' * R_data_arr( :, :, i );
end


R_arr_save = zeros( 3, 3*N );
R_arr_del_save = zeros( 3, 3*N );

for i = 1 : N
    R_arr_save( :, 3*(i-1)+1:3*i )     = R_data_arr( :, :, i );
    R_arr_del_save( :, 3*(i-1)+1:3*i ) = R_arr_del( :, :, i );
end

% At the end-of the day, what is matters is the delta movement
% Saving the orientation data as an SO(3) matrix
writematrix(  R_arr_save    , [ dir_name, 'orientation_data.csv'     ] ) 
writematrix(  R_arr_del_save, [ dir_name, 'orientation_del_data.csv' ] ) 


