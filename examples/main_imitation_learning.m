%% [Example Script] Imitation Learning.
% [Author] Moses Chong-ook Nah
% [Date]   2023.08.15

%% [Initialization] 
clear; close all; clc;
fig_config( 'fontSize', 20, 'markerSize', 10 )

%% [Imitation Learning]
% Requires q, dq, ddq for imitation learning

% Minimum jerk Trajectory
% Sample data for q, dq, ddq
q0i = 1.0;
q0f = 2.0;
D   = 2.0;
t0i = 0.5;
N   = 40;

% Parameters of the Transformation System
alpha_z = 1.0;
alpha_s = 1.0;
beta_z  = alpha_z/4;
tau     = D;
g       = q0f;
y0      = q0i;
z0      = 0.0;

trans_sys = TransformationSystem( alpha_z, beta_z, D, y0, z0 );

% Canonical System
cs = CanonicalSystem( 'discrete', tau, alpha_s );
fs = NonlinearForcingTerm( cs, N );

%% Conduct Imitation Learning

% Assume that we have 20 sample points for the minimum jerk trajectory.
P = 20;

t_P = linspace( 0.0, D, P );
  y_des_arr = zeros( 1, P );
 dy_des_arr = zeros( 1, P );
ddy_des_arr = zeros( 1, P );

% Learning the N weights based on Locally Weighted Regression
w_arr = zeros( 1, N );

for i = 1: N
    a_arr   = zeros( 1, P );
    f_arr   = zeros( 1, P );
    phi_mat = zeros( P, P );
    
    for j = 1 : P
        a_arr( j ) = ( g - y0 ) * cs.calc( t_P( j ) );
        
        [ y_des, dy_des, ddy_des ] = min_jerk_traj( t_P( j ), q0i, q0f, D, 0 );
          y_des_arr( j ) =   y_des;
         dy_des_arr( j ) =  dy_des;
        ddy_des_arr( j ) = ddy_des;
        
        f_arr( j ) = tau^2 * ddy_des + alpha_z * tau * dy_des + alpha_z * beta_z * ( y_des - g );
        phi_mat( j, j ) = fs.calc_ith( t_P( j ), i );
    end
    
    w_arr( i ) = ( a_arr * phi_mat * f_arr' ) / ( a_arr * phi_mat * a_arr' );
end

%% Generating the desired trajectory from the learned weight

dt = 1e-3;
Nt = 7000;

t_arr = dt * (0:Nt);

y_arr = zeros( 1, Nt + 1 );
z_arr = zeros( 1, Nt + 1 );

y_arr( 1 ) = y0;
z_arr( 1 ) = z0;

t = 0;

tmp_arr = zeros( 1, Nt );

for i = 1 : Nt
    
    % Calculating the input from the weights
    phi_act_arr = zeros( 1, N );
    for j = 1 : N
        phi_act_arr( j ) = fs.calc_ith( t-t0i, j );
    end
    
    if sum( phi_act_arr ) <= 1e-9
        f_input = 0;
    else
        f_input = sum( phi_act_arr .* w_arr )/sum( phi_act_arr ) * cs.calc( t-t0i ) * ( g - y0 );
    end
    
    tmp_arr( i ) = sum( phi_act_arr );
    
    [ y, z, ~, ~ ] = trans_sys.step( g, f_input, dt );
    y_arr( i + 1 ) = y;
    z_arr( i + 1 ) = z;
    
    t = t + dt;
end

%% 
plot( y_arr )