% [Title]     Basic Imitation Learning
% [Author]    Moses Chong-ook Nah
% [Email]     mosesnah@mit.edu
% [Update]    At 2023.11.12

%% (--) Initialization
clear; close all; clc;
fig_config( 'fontSize', 20, 'markerSize', 10 )

%% =======================================================
%% (1-) Imitation Learning for Minimum-jerk Trajectory
%%  -- (1A) Initialization

% The trajectory we aim to imitate is the Minimum jerk Trajectory
% The parameters of the trajectory include:
%   (1) Initial position (y0d) (column vec.)
%   (2)   Final position (gd) (column vec.)
%   (3) Duration (D) 
%   (4) Starting time (t0i)
% For (1), (2) the g subscript is to emphasize the initial/final
% position of a "d"emonstrated trajectory

y0d = [ 4.0; 4.0 ];
gd  = [ 0.0; 0.0 ];
D   = 1.0;
n   = length( y0d );

% Parameters for the Canonical System
tau     =    D;
alpha_s =  1.0;

% The number of Basis Function for the Nonlinear Forcing Term
N = 50;

% Parameters of the Transformation System
alpha_z   = 10.0;
beta_z    = 0.25 * alpha_z;

cs        = CanonicalSystem( 'discrete', tau, alpha_s );
trans_sys = TransformationSystem( alpha_z, beta_z, cs );
fs        = NonlinearForcingTerm( cs, N );

%%  -- (1B) Learning the weights via Linear Least Square Regression

% P sample points for the minimum jerk trajectory.
P = 100;

% Equal Sampling along time with duration D
t_P = linspace( 0.0, D, P );


[ y_des, dy_des, ddy_des ] = min_jerk_traj( y0d, gd, D, t_P, 0 );


% Calculating the required Nonlinear Forcing Term
f_arr   = trans_sys.get_desired( y_des, dy_des, ddy_des, gd );

% The phi matrix 
Phi_mat = zeros( P, N );

% Interating along the sample points
for i = 1 : P 
    t = t_P( i );
    Phi_mat( i, : ) = fs.calc_multiple_ith( t, 1:N )/ fs.calc_whole_at_t( t ) * cs.calc( t );
end

% Scaling matrix
% There are two options, from Koutras and Doulgeri (2020) and the Ijspeert et al. (2013).
% For this time, we use Ijspeert et al. (2013).
scl = diag( gd - y0d );
w_arr = transpose( ( Phi_mat' * Phi_mat )^(-1) * Phi_mat' * f_arr' * scl^(-1) );

%%  -- (1C) Rollout Generating a Full trajectory with the Transformation System

% Rollout with the weight array 
t0i   = 0.0;
T     = 8.0;
Nt    = 10000;
t_arr = linspace( 0, T, Nt+1 ); 

% The new y0, z0, g
y0 = y0d;
z0 = zeros( n, 1 );
g  = zeros( n, 1 );

input_arr_discrete = fs.calc_forcing_term( t_arr( 1:end-1 ), w_arr, t0i, diag( g - y0 ) );

[ y_arr, ~, dy_arr ] = trans_sys.rollout( y0, z0, g, input_arr_discrete, t0i, t_arr  );

f = figure( ); a = axes( 'parent', f );
hold( a, 'on' )
plot( a, t_arr, y_arr' , 'linewidth', 3 )
plot( a, t_arr, dy_arr', 'linewidth', 3 )

%% -- (1D) Saving the Weights and Parameters of the Learned Trajectory 

% All the necessary data 
data           = struct;
data.alpha_z   = alpha_z;
data.beta_z    = beta_z;
data.cs        = cs;
data.trans_sys = trans_sys;
data.fs        = fs;
data.g         = g;
data.w_arr     = w_arr;

% save( './learned_parameters/min_jerk_traj' , 'data' );

%% (2-) Imitation Learning for Actual Learning Trajectory
%%  --  (2A) Definition of Trajectories

t_sym = sym( 't_sym' );

l  = 2;
D  = 3.0;
v  = l/D;        
tn = t_sym/D;        

% Position
px = 0;
py =   l * ( 10 * tn^3 - 15 * tn^4 + 6 * tn^5 );
pz = 0.5 * ( cos( 2*2*pi/l * ( l * tn ) ) - 1 );

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
N = 100;
cs        = CanonicalSystem( 'discrete', D, alpha_s );
trans_sys = TransformationSystem( alpha_z, beta_z, cs );
fs        = NonlinearForcingTerm( cs, N );

%%  --  (2B) Learning the weights of the trajectory with LLS

% Calculate the f_array 
f_arr = trans_sys.get_desired( y_des_arr, dy_des_arr, ddy_des_arr, g_d ); 

% Learning the weights with least-square method 
n = length( g_d );

% A More simple calculation of the weights 
Phi_mat = zeros( P, N );

for i = 1 : P 
    Phi_mat( i, : ) = fs.calc_ith_arr( t_P( i ), 1:N )/ fs.calc_whole_at_t( t_P( i ) ) * cs.calc( t_P( i ) );
end

tmp = 1./( g - y0 );
tmp( isinf( tmp ) | isnan( tmp ) ) = 0; 

w_arr_LSS = transpose( (Phi_mat' * Phi_mat)^(-1) * Phi_mat' * f_arr' * diag( tmp ) );


for k = 1 : n
    phi_mat = zeros( P, N );
    
    for i = 1 : P
        for j = 1 : N
            phi_mat( i, j ) = fs.calc_ith( t_P( i ), j ) / fs.calc_whole_at_t( t_P( i ) ) * cs.calc( t_P( i ) );
        end
    end
    
    % Get w_arr with Least square solution
    w_arr_LSS( k, : ) = transpose( ( phi_mat' * phi_mat )^(-1) * phi_mat' * f_arr( k, : )' );

end

%%  --  (2C) Rollout with the weighting matrices

% Rollout with the weight array 
t0i = 1.0;
T   = 5;
N   = 3000;
t_arr = linspace( 0, T, N+1 ); 

y0_new = y0_d;
z0_new = zeros( n, 1 );

f = figure( ); a = axes( 'parent', f );
hold on

rot1 = roty( 30 );
tmp_scl = 1.0;

for angle = 0:60:360
    g_new = tmp_scl * rotz( angle ) * rot1* g_d;
    input_arr = fs.calc_forcing_term( t_arr( 1:end-1 ), w_arr_LSS, t0i, tmp_scl * rotz( angle ) * rot1 );
    [ y_arr, z_arr, dy_arr ] = trans_sys.rollout( y0_new, z0_new, g_new, input_arr, t0i, t_arr  );    
    plot3( a, y_arr( 1, : ), y_arr( 2, : ), y_arr( 3, : ), 'linewidth', 3, 'color', 'k' )
end

tmp_scl = 2.0;
rot1 = roty( 60 );
y0_new = [ 1.0, 1.0, 1.0 ]';
for angle = 0:60:360
    g_new = tmp_scl * rotz( angle ) * rot1* g_d;
    input_arr = fs.calc_forcing_term( t_arr( 1:end-1 ), w_arr_LSS, t0i, tmp_scl * rotz( angle ) * rot1 );
    [ y_arr, z_arr, dy_arr ] = trans_sys.rollout( y0_new, z0_new, g_new, input_arr, t0i, t_arr  );    
    plot3( a, y_arr( 1, : ), y_arr( 2, : ), y_arr( 3, : ), 'linewidth', 3, 'color', 'b' )
end

axis equal
lw = 2;
set( a, 'view',  [38.8448, 16.6583], 'xlim', [-lw, lw], 'ylim', [-lw, lw], 'zlim', [-lw, lw] )

%% =======================================================
%% (3-) Imitation Learning for Oscillation
%%  -- (1A) Oscillation

% The number of Basis Function for the Nonlinear Forcing Term
N = 20;

% Parameters for the Canonical System
T_p     = 2.0;
tau     = T_p/(2*pi);
alpha_s =  1.0;

% Parameters of the Transformation System
alpha_z   = 100.0;
beta_z    = 0.5 * alpha_z;

cs        = CanonicalSystem( 'rhythmic', tau, alpha_s );
trans_sys = TransformationSystem( alpha_z, beta_z, cs );
fs        = NonlinearForcingTerm( cs, N );

% P sample points for the minimum jerk trajectory.
P = 200;

% Equal Sampling along time with duration D
t_P = linspace( 0.0, T_p, P );

% Desired Trajectories, position, velocity and acceleration.
  y_des = zeros( 2, P );
 dy_des = zeros( 2, P );
ddy_des = zeros( 2, P );

  y_des( 1, : ) = cos( pi * t_P );
  y_des( 2, : ) = sin( pi * t_P );

 dy_des( 1, : ) = -pi * sin( pi * t_P );
 dy_des( 2, : ) =  pi * cos( pi * t_P );

ddy_des( 1, : ) = -pi^2 * cos( pi * t_P );
ddy_des( 2, : ) = -pi^2 * sin( pi * t_P );

% Calculating the required Nonlinear Forcing Term
f_arr   = trans_sys.get_desired( y_des, dy_des, ddy_des, zeros( 2, 1) );

% The phi matrix 
Phi_mat = zeros( P, N );

% Interating along the sample points
for i = 1 : P 
    t = t_P( i );
    Phi_mat( i, : ) = fs.calc_multiple_ith( t, 1:N )/ fs.calc_whole_at_t( t );
end

% Scaling matrix
% There are two options, from Koutras and Doulgeri (2020) and the Ijspeert et al. (2013).
% For this time, we use Ijspeert et al. (2013).
w_arr = transpose( ( Phi_mat' * Phi_mat )^(-1) * Phi_mat' * f_arr' );

% Rollout with the weight array 
t0i   = 0.0;
T     = 4*T_p;
Nt    = 10000;
t_arr = linspace( 0, T, Nt+1 ); 

% The new y0, z0, g
y0 =  y_des( :, 1 );
z0 = dy_des( :, 1 )*tau;
g  = zeros( 2, 1 );

input_arr_rhythmic = fs.calc_forcing_term( t_arr( 1:end-1 ), w_arr, t0i, eye( 2 ) );

[ y_arr, z_arr, dy_arr ] = trans_sys.rollout( y0, z0, g, input_arr_rhythmic, t0i, t_arr  );

f = figure( ); a = axes( 'parent', f );
hold( a, 'on' )
plot( a, y_arr( 1, : ), y_arr( 2, : ) )
axis equal

%%  -- (2B) Mixing the Two Movement Primitives

alpha = 0:0.1:1.0;
f = figure( ); a = axes( 'parent', f );
hold on, axis equal
for i = 1:length( alpha )
    a = alpha( i );
    [ y_arr, z_arr, dy_arr ] = trans_sys.rollout( a*y0+(1-a)*y0d, a*z0, g, a * input_arr_rhythmic + (1-a)*input_arr_discrete, t0i, t_arr  );
    plot( 2*i + y_arr( 1, : ), y_arr( 2, : ) )
end