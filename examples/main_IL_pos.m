%% [Example Script] Imitation Learning in R3 position Space
% [Author] Moses Chong-ook Nah
% [Date]   2023.08.15

%% (--) [Initialization] 
clear; close all; clc;
fig_config( 'fontSize', 20, 'markerSize', 10 )

%% (1-) [Imitation Learning]
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
g       = q0f;
y0      = q0i;
z0      = 0;

cs        = CanonicalSystem( 'discrete', tau, alpha_s );
trans_sys = TransformationSystem( alpha_z, beta_z, tau, y0, z0 );
fs        = NonlinearForcingTerm( cs, N );

%% ---- (1B) Learning weights via Locally Weighted Regression
is_plot = false;

% Assume that we have P sample points for the minimum jerk trajectory.
P = 100;

t_P = linspace( 0.0, D, P );
  y_des_arr = zeros( 1, P );
 dy_des_arr = zeros( 1, P );
ddy_des_arr = zeros( 1, P );

% Learning the N weights based on Locally Weighted Regression
w_arr = zeros( 1, N );

for i = 1: N
    a_arr   = zeros( 1, P );
    f_arr   = zeros( 1, P );
    phi_arr = zeros( 1, P );
    
    for j = 1 : P
        a_arr( j ) = ( g - y0 ) * cs.calc( t_P( j ) );
        
        [ y_des, dy_des, ddy_des ] = min_jerk_traj( t_P( j ), q0i, q0f, D, 0 );
          y_des_arr( j ) =   y_des;
         dy_des_arr( j ) =  dy_des;
        ddy_des_arr( j ) = ddy_des;
        
        f_arr( j ) = tau^2 * ddy_des + alpha_z * tau * dy_des + alpha_z * beta_z * ( y_des - g );
        phi_arr( j ) = fs.calc_ith( t_P( j ), i );
    end
    
    w_arr( i ) = sum( a_arr .* phi_arr .* f_arr ) / sum( a_arr .* phi_arr .* a_arr );
end

if is_plot
   
    subplot( 2, 1, 1 )
    hold on
    % Plotting the y_des, dy_des, ddy_des of the imitated trajectory
    plot( t_P,   y_des_arr );
    plot( t_P,  dy_des_arr );
    plot( t_P, ddy_des_arr );
    
    subplot( 2, 1, 2 );
    hold on
    for i = 1 : N
        plot( t_P, w_arr( i ) * fs.calc_ith( cs.calc( t_P ), i ) )
    end
end

%% ---- (1C) Generating a Full trajectory with the Transformation System

% The time step of the simulation and its number of iteration
dt = 1e-3;
Nt = 2000;

% The total time and its time array
T  = dt * Nt;
t_arr = dt * (0:Nt);

% For plotting and saving the data
y_arr = zeros( 1, Nt + 1 );
z_arr = zeros( 1, Nt + 1 );

y_arr( 1 ) = y0;
z_arr( 1 ) = z0;

t = 0;


f_input_arr = zeros( 1, Nt );
for i = 1 : Nt
    
    % Before conducting the movement
    % Maintaining that posture
    if t <= t0i
        y_arr( i + 1 ) = y0;
        z_arr( i + 1 ) = z0;
        
    % During the movement
    elseif t0i <= t
        
        if t<= t0i + D

            % taking off the initial time offset
            t_tmp = t - t0i

            % Calculating the input from the weights
            phi_act_arr = zeros( 1, N );
            for j = 1 : N
                phi_act_arr( j ) = fs.calc_ith( t_tmp, j );
            end

            if sum( phi_act_arr ) == 0
                f_input = 0;
            else
                f_input = sum( phi_act_arr .* w_arr )/sum( phi_act_arr ) * cs.calc( t_tmp ) * ( g - y0 );
            end
            
        else
            f_input = 0; 
        end

        f_input_arr( i ) = f_input; 
        [ y, z, ~, ~ ] = trans_sys.step( g, f_input, dt );
        y_arr( i + 1 ) = y;
        z_arr( i + 1 ) = z;
    end

    t = t + dt;
end

%% Plotting the Data

is_plot = true;

if is_plot 
    hold on
    plot( t_arr, y_arr )
    set( gca, 'xlim', [0, T] )
end
