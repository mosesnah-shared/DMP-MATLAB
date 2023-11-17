% [Title]     Data Filter and learn the Weight Array, Orientation, Unit Quaternion
% [Author]    Moses Chong-ook Nah
% [Email]     mosesnah@mit.edu
% [Update]    At 2023.11.12

% [Notes]
% Method heavily based on the following paper by Koutras and Doulgeri (2020)
% K.Leonidas and Z.Doulgeri. "A correct formulation for the orientation dynamic movement primitives for robot control in the cartesian space."

%% (--) Initialization
clear; close all; clc;
fig_config( 'fontSize', 20, 'markerSize', 10 )

%% =======================================================
%% (1-) Orientation Data Filtering and Learning the Weights
%%  -- (1A) Calling the data and get the Forward Kinematics 

% All dataset are saved under 'data' directory
idx = 2;
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
% We will learn the orientation, but also saving the position for backup
p_arr_raw    = zeros( 3, Np );

% First and foremost, we filter out the q_arr and dq_arr 
q_arr_filt    = zeros( 7, Np );
dq_arr_filt   = zeros( 7, Np );


% Filtered joint-position and velocity
for i = 1 : 7
    % Gaussian Filtering of Position 
     q_arr_filt( i, : ) = smoothdata( q_arr_demo( i, : ), "gaussian", 50 );

    % Diff and Gaussian Filter of Velocity  
    dq_arr_filt( i, : ) = data_diff( q_arr_filt( i, : ) );
    dq_arr_filt( i, : ) = smoothdata( dq_arr_filt( i, : ), "gaussian", 50 );

end


% Calculating the position/orientation with Forward Kinematics
% Using the filtered q and dq values
p_arr_filt    = zeros( 3, Np );
quat_arr_filt = zeros( 4, Np );
R_arr_filt    = zeros( 3, 3, Np );

for i = 1 : Np
    H = robot.getForwardKinematics( q_arr_filt( :, i ) );
    p_arr_filt( :, i ) = H( 1:3, 4   );
    R_arr_filt( :, :, i ) = H( 1:3, 1:3 ); 
    quat_arr_filt( :, i ) = rotm2quat( H( 1:3, 1:3 ) );
end

% [Note] [Moses C. Nah] [2023.11.13]
% For changing the Rotation Matrix to Quaternion
% The great references for this are shown below:
% [1] Appendix B2, Allmendinger, Felix. "Computational methods for the kinematic analysis of diarthrodial joints" (2015)
% [2] [1] is from Shepperd, Stanley W. "Quaternion from rotation matrix." (1978):

% Double check the results 
f = figure( ); a = axes( 'parent', f );
hold( a, 'on' ); axis equal
set( a, 'view', [58.6300, 5.5740] )
plot3( a, p_arr_filt( 1, : ), p_arr_filt( 2, : ), p_arr_filt( 3, : ), 'linewidth', 3, 'color', 'black' )
scatter3( a, p_arr_filt( 1,   1 ), p_arr_filt( 2,   1 ), p_arr_filt( 3,   1 ), 500, 'filled', 'o', 'markeredgecolor', [0.0000, 0.4470, 0.7410], 'markerfacecolor', 'white', 'linewidth', 5 )
scatter3( a, p_arr_filt( 1, end ), p_arr_filt( 2, end ), p_arr_filt( 3, end ), 500, 'filled', 'o', 'markeredgecolor', [0.8500, 0.3250, 0.0980], 'markerfacecolor', 'white', 'linewidth', 5 )

tmp_scl = 0.02;
for i = 1 : 100: Np
    
    scatter3( a, p_arr_filt( 1, i ), p_arr_filt( 2, i ), p_arr_filt( 3, i ), 500, 'filled', 'o', 'markeredgecolor', 'black', 'markerfacecolor', 'white', 'linewidth', 5 )    
     quiver3( a, p_arr_filt( 1, i ), p_arr_filt( 2, i ), p_arr_filt( 3, i ), tmp_scl * R_arr_filt( 1, 1, i ), tmp_scl * R_arr_filt( 2, 1, i ), tmp_scl * R_arr_filt( 3, 1, i ), 'linewidth', 4, 'color', 'r' )
     quiver3( a, p_arr_filt( 1, i ), p_arr_filt( 2, i ), p_arr_filt( 3, i ), tmp_scl * R_arr_filt( 1, 2, i ), tmp_scl * R_arr_filt( 2, 2, i ), tmp_scl * R_arr_filt( 3, 2, i ), 'linewidth', 4, 'color', 'g' )
     quiver3( a, p_arr_filt( 1, i ), p_arr_filt( 2, i ), p_arr_filt( 3, i ), tmp_scl * R_arr_filt( 1, 3, i ), tmp_scl * R_arr_filt( 2, 3, i ), tmp_scl * R_arr_filt( 3, 3, i ), 'linewidth', 4, 'color', 'b' )

end

%%  -- (1B) Pre-processing the Quaternion Data
% What we need for learning is the error quaternion
% In detail, given the goal quaternion qg, the error quaternion is 
% qe = 2Im( Log( conj( q(t) ) * qg ) );
% Note that we have switched the location of qg

% Calculate the error array, given the quaternion goal
quat_goal = quat_arr_filt( :, end );
error_arr = zeros( 3, Np );

for i = 1 : Np
    error_arr( :, i ) = get_quat_error( quat_arr_filt( :, i ), quat_goal  );
end

% For training, we need to calculate the error array value 
% First, we need to calculate the time derivative of quaternion
% For that, we need the spatial angular velocity 

ws_arr    = zeros( 3, Np );
dquat_arr = zeros( 4, Np );
for i = 1 : Np
    JS  = robot.getSpatialJacobian( q_arr_filt( :, i ) );
    JSr = JS( 4:6, : );
    ws_arr( :, i ) = JSr * dq_arr_filt( :, i );

    dquat_arr( :, i ) = 1/2 * quat_mul( R3_to_quat( ws_arr( :, i ) ), quat_arr_filt( :, i  )' );
end

dLogquat_arr = zeros( 3, Np );

for i = 1 : Np
    dLogquat_arr( :, i ) = dLogQuat( quat_arr_filt( :, i ), dquat_arr( :, i ) );
end

% Finally, calculate the ddLog part by simply time derivative
ddLogquat_arr = zeros( 3, Np );

% Filtered joint-position and velocity
for i = 1 : 3

    % Diff and Gaussian Filter of Velocity  
    ddLogquat_arr( i, : ) = data_diff( dLogquat_arr( i, : ) );
    ddLogquat_arr( i, : ) = smoothdata( ddLogquat_arr( i, : ), "gaussian", 50 );

end

%%  -- (1C) Learning the weights by LLS regression
% The form is exactly the same with Linear DMP
% Parameters of DMP
alpha_s   =  2.0;
alpha_z   = 2000.0;
beta_z    = 0.5 * alpha_z;
N         = 100;
tau       = max( t_arr_demo );

% The Three elements of DMP
cs        = CanonicalSystem( 'discrete', tau, alpha_s );
trans_sys = TransformationSystemQuat( alpha_z, beta_z, cs );
fs        = NonlinearForcingTerm( cs, N );

f_arr = trans_sys.get_desired( error_arr, dLogquat_arr, ddLogquat_arr );

% The phi matrix 
P = length( t_arr_demo );
Phi_mat = zeros( P, N );

% Calculate the Phi matrix for Weight Learning
% Note that compared to Ijspeert 2013, this followed Koutras and Doulgeri (2020)
for i = 1 : P 
    t = t_arr_demo( i );
    Phi_mat( i, : ) = fs.calc_multiple_ith( t, 1:N )/ fs.calc_whole_at_t( t ) * cs.calc( t );
end

% Scaling Matrix
Sr = diag( get_quat_error( quat_arr_filt( :, 1 ), quat_arr_filt( :, end ) ) );

% Learning the weights with Linear Least-square fitting
w_arr = transpose( ( Phi_mat' * Phi_mat )^(-1) * Phi_mat' * f_arr' * inv( Sr ) );

% Double check the result
% Need the rollout
% Double check the Learned weights with basic DMP
% Initial conditions
eq0 =    error_arr( :, 1 );
z0  = dLogquat_arr( :, 1 )/tau;
t0i = 0.0;
T   = max( t_arr_demo );
dt  = 1e-4;
t_arr = 0:dt:T;
Ns = length( t_arr );

% Calculate the nonlinear forcing term
input_arr = fs.calc_forcing_term( t_arr( 1:end-1 ), w_arr, t0i, Sr  );

[ eq_arr, z_arr, deq_arr ] = trans_sys.rollout( eq0, z0, 1*input_arr, t0i, t_arr  );

f = figure( ); a = axes( 'parent', f );
hold( a, 'on' )

for i = 1: 3
    plot( t_arr_demo, error_arr( i, : ), 'linestyle', '--', 'color', 'k', 'linewidth', 10 )
    plot( t_arr     ,    eq_arr( i, : ), 'linestyle', '-' , 'color', 'r', 'linewidth', 5  )
end

%%  -- (1D) Revive the Quaternion from the error quaternion

quat_arr_gen = zeros( 4, Ns );
R_arr_gen    = zeros( 3, 3, Ns );
R_arr_cpr    = zeros( 3, 3, Ns );
quat_g = quat_arr_filt( :, end );

for i = 1 : Ns
    tmp = R3_to_quat( eq_arr( :,i ) );
    quat_arr_gen( :, i ) = quat_mul( quat_g', quat_conj( ExpQuat( 0.5 * tmp ) ) );

    R_arr_gen( :, :, i ) = quat_to_SO3( quat_arr_gen( :, i ) );
    
    ttmp = quat_arr_gen( :, i );
    R_arr_cpr( :, :, i ) = quat2rotm( quaternion( ttmp( 1 ), ttmp( 2 ), ttmp( 3), ttmp( 4 ) ) );
end


% Double check the results 
f = figure( ); a = axes( 'parent', f );
hold( a, 'on' ); axis equal
set( a, 'view', [58.6300, 5.5740] )
plot3( a, p_arr_filt( 1, : ), p_arr_filt( 2, : ), p_arr_filt( 3, : ), 'linewidth', 3, 'color', 'black' )
scatter3( a, p_arr_filt( 1,   1 ), p_arr_filt( 2,   1 ), p_arr_filt( 3,   1 ), 500, 'filled', 'o', 'markeredgecolor', [0.0000, 0.4470, 0.7410], 'markerfacecolor', 'white', 'linewidth', 5 )
scatter3( a, p_arr_filt( 1, end ), p_arr_filt( 2, end ), p_arr_filt( 3, end ), 500, 'filled', 'o', 'markeredgecolor', [0.8500, 0.3250, 0.0980], 'markerfacecolor', 'white', 'linewidth', 5 )

tmp_scl = 0.02;
for i = 1 : 100: Np

    
    scatter3( a, p_arr_filt( 1, i ), p_arr_filt( 2, i ), p_arr_filt( 3, i ), 500, 'filled', 'o', 'markeredgecolor', 'black', 'markerfacecolor', 'white', 'linewidth', 5 )    
     quiver3( a, p_arr_filt( 1, i ), p_arr_filt( 2, i ), p_arr_filt( 3, i ), tmp_scl * R_arr_filt( 1, 1, i ), tmp_scl * R_arr_filt( 2, 1, i ), tmp_scl * R_arr_filt( 3, 1, i ), 'linewidth', 4, 'color', 'r' )
     quiver3( a, p_arr_filt( 1, i ), p_arr_filt( 2, i ), p_arr_filt( 3, i ), tmp_scl * R_arr_filt( 1, 2, i ), tmp_scl * R_arr_filt( 2, 2, i ), tmp_scl * R_arr_filt( 3, 2, i ), 'linewidth', 4, 'color', 'g' )
     quiver3( a, p_arr_filt( 1, i ), p_arr_filt( 2, i ), p_arr_filt( 3, i ), tmp_scl * R_arr_filt( 1, 3, i ), tmp_scl * R_arr_filt( 2, 3, i ), tmp_scl * R_arr_filt( 3, 3, i ), 'linewidth', 4, 'color', 'b' )

     % Find the time for this array
     tt = t_arr_demo( i );

     % Find the closest index
     [ d, ix ] = min( abs( t_arr-tt ) );


     quiver3( a, p_arr_filt( 1, i ), p_arr_filt( 2, i ), p_arr_filt( 3, i ), tmp_scl * R_arr_cpr( 1, 1, ix ), tmp_scl * R_arr_cpr( 2, 1, ix ), tmp_scl * R_arr_cpr( 3, 1, ix ), 'linewidth', 4, 'color', 'r' )
     quiver3( a, p_arr_filt( 1, i ), p_arr_filt( 2, i ), p_arr_filt( 3, i ), tmp_scl * R_arr_cpr( 1, 2, ix ), tmp_scl * R_arr_cpr( 2, 2, ix ), tmp_scl * R_arr_cpr( 3, 2, ix ), 'linewidth', 4, 'color', 'g' )
     quiver3( a, p_arr_filt( 1, i ), p_arr_filt( 2, i ), p_arr_filt( 3, i ), tmp_scl * R_arr_cpr( 1, 3, ix ), tmp_scl * R_arr_cpr( 2, 3, ix ), tmp_scl * R_arr_cpr( 3, 3, ix ), 'linewidth', 4, 'color', 'b' )

end

%%  -- (1E) Saving the Data, the weights for the learned trajectory

data = struct;

% Parameters of DMP and Learned Weighted
data.name    = 'draw_M';
data.tau     = tau;
data.alpha_s = alpha_s;
data.alpha_z = alpha_z;
data.beta_z  =  beta_z;
data.weight  =   w_arr;

% Initial and Goal Location
data.goal = quat_arr_filt( :, end );
data.z0   = dLogquat_arr( :, 1 )/tau;

save( [ 'learned_parameters/', traj_names{ idx } , '_orient.mat' ], 'data' );

%% =======================================================
%% (2-) Extension
%%  -- (2A) Visualizing the Learned Rotation

% Create a star-like object
% [REF] https://www.mathworks.com/help/matlab/ref/hgtransform.html
f = figure(  ); a = axes( 'parent', f );
axis equal; 
[ x, y, z] = cylinder( [0.2 0.0] );
h( 1 ) = surface( a,  x,  y,  z, 'FaceColor', 'red'     );
h( 2 ) = surface( a,  x,  y, -z, 'FaceColor', 'green'   );
h( 3 ) = surface( a,  z,  x,  y, 'FaceColor', 'blue'    );
h( 4 ) = surface( a, -z,  x,  y, 'FaceColor', 'cyan'    );
h( 5 ) = surface( a,  y,  z,  x, 'FaceColor', 'magenta' );
h( 6 ) = surface( a,  y, -z,  x, 'FaceColor', 'yellow'  );
set( a, 'xlim', [ -2.0, 2.0], 'ylim', [ -2.0, 2.0 ], 'zlim', [ -2.0, 2.0 ] )
view( a, 3 )

t = hgtransform( 'Parent' ,a );
set( h, 'Parent', t )

% Saving the video for the rotations
v = VideoWriter( 'rotation_video.mp4' );
open(v)

Nfs = round( 1/dt );
Nstep = round( Nfs/30 );
for i = 1 : Nstep : Ns
    tmp = eye( 4 );
    tmp( 1:3, 1:3 ) = R_arr_gen( :, :, i );
    set( t ,'Matrix', tmp )
    drawnow 
    
    frame = getframe(f);
    writeVideo( v, frame )

end
close( v )


