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
N  = 100;

% Parameters of the 3 DMPs
alpha_z = 160.0;
alpha_s = 20.0;
beta_z  = 10.0;
tau     = max( t_arr );


% Get the number of total samples
[ ~, P_total ] = size( quat );

idx = 1 : 5 : P_total;
P = length( idx );

% Get the goal location  
g     = quat( :, idx( end ) );
cs    = CanonicalSystem( 'discrete', tau, alpha_s );
fs_wx = NonlinearForcingTerm( cs, N );
fs_wy = NonlinearForcingTerm( cs, N );
fs_wz = NonlinearForcingTerm( cs, N );

fs_w = cell( 1, 3 ); 
fs_w{ 1 } = fs_wx; fs_w{ 2 } = fs_wy; fs_w{ 3 } = fs_wz;

% Learning the N weights based on Locally Weighted Regression
w_arr = zeros( 3, N );

tmp0 = quat_mul( g, quat_conj( quat( :, 1 ) ) ); 
vec0 = tmp0( 2:4 );

for k = 1 : 3

    for i = 1: N

        a_arr   = zeros( 1, P );
        f_arr   = zeros( 1, P );
        phi_arr = zeros( 1, P );

        for j = 1 : P

            tmp1 = quat_mul( g, quat_conj( quat( :, idx( j ) ) ) ); 
            vec1 = tmp1( 2 : 4 );

            ddy_des = domega_arr( :, idx( j ) );
             dy_des =  omega_arr( :, idx( j ) );

            a_arr( j ) = vec0( k ) * cs.calc( t_arr( idx( j ) ) );
            f_arr( j ) = tau * ddy_des( k ) + beta_z * dy_des( k ) - alpha_z * vec1( k );
            phi_arr( j ) = fs_w{ k }.calc_ith( t_arr( idx( j ) ), i );
            
        end

        if( sum( a_arr .* phi_arr .* a_arr ) ~= 0 )
            w_arr( k, i ) = sum( a_arr .* phi_arr .* f_arr ) / sum( a_arr .* phi_arr .* a_arr );
        else
            w_arr( k, i ) = 0;
        end

    end

end

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
    
    tmp1 = quat_mul( g, quat_conj( y ) ); 
    vectmp = tmp1( 2: 4 );
    
    for k = 1 : 3
        
        phi_act_arr = zeros( 1, N );
        
        for j = 1 : N
            phi_act_arr( j ) = fs_w{ k }.calc_ith( t, j );
        end

        if sum( phi_act_arr ) == 0
            f_input( k ) = 0;
        else
            f_input( k ) = sum( phi_act_arr .* w_arr( k, : ) )/sum( phi_act_arr ) * cs.calc( t ) * vectmp( k );
        end
    
    end

    [ quat_new, w_new, ~, ~ ] = trans.step( g, zeros( 3, 1), dt );
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
