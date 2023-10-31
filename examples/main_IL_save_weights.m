%% ==================================================================
%% [Title] Imitation Learning and Saved the Weight Array
% Author: Moses Chong-ook Nah
%  Email: mosesnah@mit.edu
%   Date: 2023.10.18
%% ==================================================================

%% [0A] Initialization
clear; close all; clc;
fig_config( 'fontSize', 20, 'markerSize', 10 )

%% [1-] Imitation Learning of a Given Trajectory
%% ---- [1A] Define Trajectory

is_check = true; 

% Define the Analytical Trajectory 2D that we aim to learn
syms t_sym

% Name of the Trajectory which we aim to learn
% There are three types:
% [1] min_jerk
% [2] cosine
% [3] radial
name_traj = 'min_jerk';

% Position
switch name_traj

    case 'min_jerk'
        D   = 3.0;
        p0i = [ 0.0; 0.00; 0.00 ];
        p0f = [ 0.0; 0.20; 0.20 ];      
        tn  = t_sym/D;

        p_func = p0i + ( p0f - p0i ) * ( 10 * tn^3 - 15 * tn^4 + 6 * tn^5 );

        px = p_func( 1 );
        py = p_func( 2 );
        pz = p_func( 3 );

    case 'cosine'
        l  = 1.0;
        D  = 3.0;
        v  = l/D;        
        tn = t_sym/D;        

        px = 0;
        py =   l * ( 10 * tn^3 - 15 * tn^4 + 6 * tn^5 );
        pz = 0.5 * ( cos( 2*2*pi/l * ( l * tn ) ) - 1 );

    case 'radial'
        r0i = 0.45;
        r0f = 0.25;
        D   = 3.0;
        tn  = t_sym/D;         

        rt    = r0i + ( r0f - r0i ) * ( -2*tn^3 + 3*tn^2 );
        theta = 2 * pi * ( -tn^4 + 2*tn^2 );

        px = 0;
        py = rt * cos( theta );
        pz = rt * sin( theta );

    otherwise
        error( 'Wrong input: %s', name_traj )
end
    
% Velocity
dpx = diff( px, t_sym );
dpy = diff( py, t_sym );
dpz = diff( pz, t_sym );

% Acceleration
ddpx = diff( dpx, t_sym );
ddpy = diff( dpy, t_sym );
ddpz = diff( dpz, t_sym );

% Saving these elements as symbolic array.
p_func   = matlabFunction( [   px;   py;   pz ] );
dp_func  = matlabFunction( [  dpx,  dpy,  dpz ] );
ddp_func = matlabFunction( [ ddpx, ddpy, ddpz ] );

% Generating the actual data from the symbolic form 
% These data will be used for Imitation Learning
dt    = 0.01;
t_arr = 0:dt:D;
P     = length( t_arr );

  p_data = zeros( 3, P );
 dp_data = zeros( 3, P );
ddp_data = zeros( 3, P );

for i = 1 : P
      p_data( :, i ) =   p_func( t_arr( i ) );
     dp_data( :, i ) =  dp_func( t_arr( i ) );
    ddp_data( :, i ) = ddp_func( t_arr( i ) );
end

yi_d = p_func( 0 ); g_d = p_func( D );

if is_check
    f = figure( ); a = axes( 'parent', f );
    plot3( a, p_data( 1, : ), p_data( 2, : ), p_data( 3, : ) )
    axis equal
end

% For Imitation Learning, one should define the number of basis function
N  = 50;

% Parameters of the 3 DMPs
% We use the identical alpha_z, beta_z values
alpha_z = 10.0;
alpha_s = 1.0;
beta_z  = 1/4 * alpha_z;
tau     = D;
z0      = 0;

cs        = CanonicalSystem( 'discrete', D, alpha_s );
trans_sys = TransformationSystem( alpha_z, beta_z, cs );
fs        = NonlinearForcingTerm( cs, N );

% Learning the weights for 3-DOF
w_arr_LSS = zeros( 3, N );


% Calculate the f_array 
f_arr = trans_sys.get_desired( p_data, dp_data, ddp_data, g_d ); 

for k = 1 : 3
    phi_mat = zeros( P, N );
    
    for i = 1 : P
        for j = 1 : N
            phi_mat( i, j ) = fs.calc_ith( t_arr( i ), j ) / fs.calc_whole_at_t( t_arr( i ) ) * cs.calc( t_arr( i ) );
        end
    end
    
    % Get w_arr with Least square solution
    w_arr_LSS( k, : ) = transpose( ( phi_mat' * phi_mat )^(-1) * phi_mat' * f_arr( k, : )' );

end

% A more simplified calculation 


%% ---- [1D] Rollout Generating a Full trajectory with the Transformation System


% Rollout with the weight array 
t0i = 1.0;
T   = 5;
N   = 3000;
t_arr = linspace( 0, T, N+1 ); 

n = length( yi_d );

y0_new = yi_d;
z0_new = zeros( n, 1 );

f = figure( ); a = axes( 'parent', f );
hold on

rot1 = roty( 30 );
tmp_scl = 1.0;

for angle = 0:120:360
    g_new  = tmp_scl * rotz( angle ) * rot1 * g_d;
    y0_new = tmp_scl * rotz( angle )* rot1 * yi_d;
    input_arr = fs.calc_forcing_term( t_arr( 1:end-1 ), w_arr_LSS, t0i, tmp_scl * rotz( angle ) * rot1 );
    [ y_arr, z_arr, dy_arr ] = trans_sys.rollout( y0_new, z0_new, g_new, input_arr, t0i, t_arr  );    
    plot3( a, y_arr( 1, : ), y_arr( 2, : ), y_arr( 3, : ), 'linewidth', 3, 'color', 'k' )
end

axis equal
lw = 2;
set( a, 'view',  [38.8448, 16.6583], 'xlim', [-lw, lw], 'ylim', [-lw, lw], 'zlim', [-lw, lw] )


%% [1F] Saving the Data

% The Weighting Matrix and additional Data
data = struct;
data.p_func   =   p_func;
data.dp_func  =  dp_func;
data.ddp_func = ddp_func;

data.name = name_traj;
data.alpha_z = alpha_z;
data.beta_z  =  beta_z;
data.alpha_s = alpha_s;

data.weight  = w_arr_LSS;
data.p0i = yi_d;
data.p0f = g_d;
 
data.cs = cs;
data.fs = fs;
data.ts = trans_sys;

save( ['learned_parameters/', name_traj ,'.mat'], 'data' );
