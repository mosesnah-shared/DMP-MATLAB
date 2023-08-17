%% [Example Script] Imitation Learning.
% [Author] Moses Chong-ook Nah
% [Date]   2023.08.15

%% [Initialization] 
clear; close all; clc;
fig_config( 'fontSize', 20, 'markerSize', 10 )
R_data = load( 'R_data.mat' );


%% [Imitation Learning]
%% ---- (1A) Parameter Initialization

% The trajectory we aim to imitate is the Minimum jerk Trajectory
% The initial (q0i), final posture (q0f), duration (D), starting time (t0i) of the trajectory
q0i = 0.0;
q0f = 1.0;
D   = 1.0;
t0i = 0.5;

% For Imitation Learning, one should define the number of basis function
N  = 20;

% Parameters of the 3 DMPs
alpha_z = 10.0;
alpha_s = 1.0;
beta_z  = 1/4 * alpha_z;
tau     = D;
g       = R_data( :, :, end );
R0      = R_data( :, :, 1   );
w0      = zeros( 1, 3 );

cs        = CanonicalSystem( 'discrete', tau, alpha_s );
trans_sys = TransformationSystem_SO3( alpha_z, beta_z, tau, R0, w0 );

%% ---- (1B) Imitation Learning

% Learning the weights of the Basis function
fs_x = NonlinearForcingTerm( cs, N );
fs_y = NonlinearForcingTerm( cs, N );
fs_z = NonlinearForcingTerm( cs, N );

% The data of R



