%% [Example Script] Imitation Learning.
% [Author] Moses Chong-ook Nah
% [Date]   2023.08.15

%% [Initialization] 
clear; close all; clc;
fig_config( 'fontSize', 20, 'markerSize', 10 )

% The SO(3) Data to learn
tmp = load( './data/R_data.mat' );
R_data = tmp.R_arr_save;
R_data = R_data * 1/0.3;
% Change the SO(3) data to Quaternions
% [REF] Allmendinger, Felix. Computational methods for the kinematic analysis of diarthrodial joints. 
%       Diss. Dissertation, Aachen, Techn. Hochsch., 2015, 2015.
[ ~, ~, N ]  = size( R_data );

quat_data = zeros( 4, N );

for i = 1 : N
    quat = SO3_to_quat( R_data( :, :, i ) ); 
    quat_data( :, i ) = quat';
end

%% [Imitation Learning]
%% ---- (1A) Parameter Initialization

% The trajectory we aim to imitate is the Minimum jerk Trajectory
% The initial (q0i), final posture (q0f), duration (D), starting time (t0i) of the trajectory
q0i = 0.0;
q0f = 1.0;
D   = 1.0;
t0i = 0.5;

% For Imitation Learning, one should define the number of basis function
N  = 50;

% Parameters of the 3 DMPs
alpha_z = 80.0;
alpha_s = 2.0;
beta_z  = 10.0;
tau     = D;
quatg   = quat_data( :, end );
quat0   = quat_data( :, 1   );
w0      = zeros( 1, 3 );

cs        = CanonicalSystem( 'discrete', tau, alpha_s );
trans_sys = TransformationSystem_quat( alpha_z, beta_z, tau, quat0, w0 );

% The time step of the simulation and its number of iteration
dt = 1e-3;
Nt = 20000;

% The total time and its time array
T  = dt * Nt;
t_arr = dt * (0:Nt);

% For plotting and saving the data
quat_arr = zeros( 4, Nt + 1 );
w_arr = zeros( 3, Nt + 1 );

quat_arr( :,1 ) = quat0;
w_arr( :,1 )    = zeros( 3, 1 ) ;

t = 0;

for i = 1 : Nt
        
    [ y, z, ~, ~ ] = trans_sys.step( quatg, zeros( 3, 1 ), dt );
    quat_arr( :,i + 1 ) = y;
    w_arr( :,i + 1 )    = z;
    
    t = t + dt;
end



tmp_color = { 'k', 'r', 'g', 'b' };

for i = 1 : 4
    plot( t_arr, quat_arr, 'color', tmp_color{ i }  )
    hold on
    plot( t_arr, repmat( quatg, size( t_arr ) ), '--', 'color', tmp_color{ i } )
end
%% ---- (1C) Plot the Quaternions

% [ ~, Ntmp ] = size( quat_arr );
% R_arr = zeros( 3, 3, Ntmp );
% 
% for i = 1 : Ntmp
%    R_arr( :, :, i ) = quat_to_SO3( quat_arr( :, i )' ); 
% end
% 
% f = figure( ); a = axes( 'parent', f );
% 
% scl = 0.3;
% scl2 = 0.00;
% hold on
% 
% for i = 1 : 20: Ntmp 
%    tmp = R_arr( :, :, i );
%    x = scl2 * i;
%    y = scl2 * i;
%    z = scl2 * i;
%    
%    quiver3( x, y, z, tmp(1,1), tmp(2,1), tmp(3,1), 'linewidth', 3, 'color', 'r' );
%    quiver3( x, y, z, tmp(1,2), tmp(2,2), tmp(3,2), 'linewidth', 3, 'color', 'g' );
%    quiver3( x, y, z, tmp(1,3), tmp(2,3), tmp(3,3), 'linewidth', 3, 'color', 'b' );
% end
