%% ==================================================================
%% [Title] Example Script - Advanced DMP from Leonidas and Doulgeri (2020)
% Author: Moses Chong-ook Nah
%  Email: mosesnah@mit.edu
%   Date: 2023.10.26
%% ==================================================================

%% [0A] Initialization
clear; close all; clc;
fig_config( 'fontSize', 20, 'markerSize', 10 )

%% [1A] Define a trajectory that we aim to learn, and also learn the weights

% Learn the trajectory residing on the yz plane
% Parameters of the trajectory
t_sym = sym( 't_sym' );

l  = 0.3;
D  = 3.0;
v  = l/D;        
tn = t_sym/D;        

% Position
px = 0;
py = l * ( 10 * tn^3 - 15 * tn^4 + 6 * tn^5 );
pz = 0.05 * ( cos( 2*2*pi/l * ( l * tn ) ) - 1 );

% Velocity
dpx = diff( px, t_sym );
dpy = diff( py, t_sym );
dpz = diff( pz, t_sym );

% Acceleration
ddpx = diff( dpx, t_sym );
ddpy = diff( dpy, t_sym );
ddpz = diff( dpz, t_sym );

% Make the trajectory as a MATLAB function
% Saving these elements as symbolic array.
p_func   = matlabFunction( [   px;   py;   pz ] );
dp_func  = matlabFunction( [  dpx,  dpy,  dpz ] );
ddp_func = matlabFunction( [ ddpx, ddpy, ddpz ] );

% Extract the sample points and learn the weights
% Assume that we have P sample points for
P = 100;

% Equal Sampling along time with duration D
t_P = linspace( 0.0, D, P );

% Desired Trajectories
  y_des_arr = zeros( 3, P );
 dy_des_arr = zeros( 3, P );
ddy_des_arr = zeros( 3, P );

for i = 1 : P
      y_des_arr( :, i ) =   p_func( t_P( i ) );
     dy_des_arr( :, i ) =  dp_func( t_P( i ) );
    ddy_des_arr( :, i ) = ddp_func( t_P( i ) );
end

% Get the initial and final position
alpha_s = 2.0;
alpha_z = 10.0;
beta_z  = 1/4 * alpha_z;
y0_d    = p_func( 0 );
g_d     = p_func( D );
z0      = dp_func( 0 )/D;

% Learn the weights using Least-square method
% Define the three elements of DMP
N = 200;
cs        = CanonicalSystem( 'discrete', D, alpha_s );
trans_sys = TransformationSystem( alpha_z, beta_z, D, y0_d, z0 );
fs        = NonlinearForcingTerm( cs, N );

% Weights arr
w_arr = zeros( 3, N );

% Learn the weights separately 
for i = 1 : 3

    f_arr   = zeros( 1, P );
    phi_mat = zeros( P, N ); 
    
    % Iterating through the sampled points
    for j = 1 : P
        f_arr( j ) = D^2 * ddy_des_arr( i, j ) + alpha_z * D * dy_des_arr( i, j ) + alpha_z * beta_z * ( y_des_arr( i, j ) - g_d( i ) );
        
        % Iterating through the basis functions
        phi_sum = fs.calc_whole_at_t( t_P( j ) );
        
        for k = 1 : N
            phi_mat( j, k ) =  fs.calc_ith( t_P( j ), k ) / phi_sum * cs.calc( t_P( j ) );
        end

    end

    w_arr( i, : ) = transpose( ( phi_mat' * phi_mat )^(-1) * phi_mat' * f_arr' );
   
end

%% [1B] Weights are learned, hence now we scale to different locations

% Generate random initial and goal position
N_traj = 5;
y0_arr = zeros( 3,  N_traj );
g0_arr = 3 * rand( 3,  N_traj );

% Also the initial condition of z0

% The starting time for each trajectory
t0_arr = ones( 1, N_traj );

dt = 1e-3;
T_arr = ( D + 2 ) * ones( 1, N_traj );
N_arr = round( T_arr/dt );

% Saving the data
y_arr  = cell( 1, N_traj );
z_arr  = cell( 1, N_traj );
dy_arr = cell( 1, N_traj );
dz_arr = cell( 1, N_traj );

% Learn the trajectories 
for i = 1 : N_traj

    % Define the y_arr for each trajectory
    y_arr{ i }  = zeros( 3, N_arr( i ) );
    z_arr{ i }  = zeros( 3, N_arr( i ) );
    dy_arr{ i } = zeros( 3, N_arr( i ) );   
    dz_arr{ i } = zeros( 3, N_arr( i ) );
    
    t = 0;
    
    % Calculate the rotation matrix and scale 
    sg = norm( g_d - y0_d )/norm( g0_arr( :, i ) - y0_arr( :, i ) );

    % Set the translational system

    % Rolling out
    for j = 1 : N_arr( i )

        % Before conducting the movement
        % Maintaining that posture
        if t <= t0_arr( i )
            y_arr{ i }( :, j + 1 ) = y0;
            z_arr{ i }( :, j + 1 ) = z0;
            
        % During the movement
        elseif t0_arr( i ) <= t
            
            if t <= t0_arr( i ) + D
    
                % taking off the initial time offset
                t_tmp = t - t0_arr( i );
    
                % Calculating the input from the weights
                % First, check if whole activation value is 0
                phi_sum = fs.calc_whole_at_t( t_tmp );
                
                f_input = 0;
    
                if phi_sum ~= 0
                    f_input = fs.calc_whole_weighted_at_t( t_tmp, w_arr_LWR )/phi_sum;

                    % Multiply f_input with the rotation matrix and etc.
                    f_input = f_input*sg*cs.calc( t_tmp );

                    % Multiply with the rotation matrix finally

                end
                
            else
                f_input = 0; 
            end
    
            [ y, z, dy, dz ] = trans_sys.step( g, f_input, dt );
            y_arr(  i + 1 ) = y;
            z_arr(  i + 1 ) = z;
            dz_arr( i + 1 ) = dz;   
            
        end
        % end of if statement for time
        
        t = t + dt;
    end
    % End of rollout 
end 
