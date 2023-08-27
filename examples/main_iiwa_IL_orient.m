%% [Example Script] Plotting the Nonlinear Forcing Terms
% [Author] Moses Chong-ook Nah
% [Date]   2023.08.15

%% [Initialization] 
clear; close all; clc;
fig_config( 'fontSize', 20, 'markerSize', 10 )

%% Load the Data that we imitate to learn 

data_raw = load( 'data/example1/iiwa_example1.mat' );

% Array of quaternion, w_arr and dw_arr 
t_arr      = data_raw.t_arr; 
quat       = data_raw.quat_arr_raw;
omega_arr  = data_raw.w_arr_filt;
domega_arr = data_raw.dw_arr_filt;

%%
% Flip the sign if detected
% For Imitation Learning, one should define the number of basis function
N  = 30;

% Parameters of the 3 DMPs
alpha_z = 8000.0;
alpha_s = 1.0;
beta_z  = 80.0;
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

    f_arr( 3*(j-1)+1:3*j ) = tau * ddy_des + beta_z * dy_des - alpha_z * quat_err;

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
Nt = P_total * 2;

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

for i = 1 : Nt
    
    % Calculating the input from the weights
    f_input = zeros( 3, 1 );
    
    vectmp = quat_log( quat_mul( g, quat_conj( y ) ) ); 
    
    % Get f_input 
    % Get the sum of phis
    phi_act_arr = zeros( 1, N );    
    for k = 1 : N
        phi_act_arr( k ) = fs_w.calc_ith( t, k );
    end
    
    f_input = zeros( 3, 1 );
    for k = 1 : N
        f_input = f_input + fs_w.calc_ith( t, k ) * w_arr( 3*(k-1)+1:3*k ) / sum( phi_act_arr ) * cs.calc( t );
    end

    [ quat_new, w_new, ~, ~ ] = trans.step( g, f_input, dt );
    quat_data_arr(  :, i  ) = quat_new;
    omega_data_arr( :, i  ) = w_new;

    y = quat_new;
    
    t = t + dt;
end


%% Change quat to SO3

[ ~, Nr] = size( quat_data_arr );

R_data_arr = zeros( 3, 3, Nr );

for i = 1 : Nr
    R_data_arr( :, :, i ) = quat_to_SO3( quat_data_arr( :, i ) );
end

tmp_color = { 'k', 'r', 'g', 'b' };
hold on
for i = 1 : 4
    plot( t_arr , quat_data_arr( i, :), 'color', tmp_color{ i } );
    plot( data_raw.t_arr, quat( i, :) , 'linestyle', '--', 'color', tmp_color{ i });
end
