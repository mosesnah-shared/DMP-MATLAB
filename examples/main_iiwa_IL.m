%% [Example Script] Plotting the Nonlinear Forcing Terms
% [Author] Moses Chong-ook Nah
% [Date]   2023.08.15

%% [Initialization] 
clear; close all; clc;
fig_config( 'fontSize', 20, 'markerSize', 10 )

%% Call iiwa from Explicit-MATLAB

is_run = false;

if is_run
    % Set figure size and attach robot to simulation
    robot = iiwa14( 'high' );
    robot.init( );

    anim = Animation( 'Dimension', 3, 'xLim', [-0.7,0.7], 'yLim', [-0.7,0.7], 'zLim', [0,1.4], 'isSaveVideo', true );
    anim.init( );
    anim.attachRobot( robot )  

    view( 90, 0 );
end

%% Load the Data that we imitate to learn 

data_raw = load( 'data/example1/iiwa_example1.mat' );

t_arr  = data_raw.t_arr;
q_arr  = data_raw.q_arr;
p_arr  = data_raw.p_arr_raw;
dp_arr = data_raw.dp_arr_filt;

[ ~, Ntmp ] = size( q_arr );


if is_run
    
    plot3(    anim.hAxes, p_arr( 1, :   ), p_arr( 2, :   ), p_arr( 3, :   ), 'linewidth', 3, 'linestyle', '--', 'color', 'k' );
    scatter3( anim.hAxes, p_arr( 1, 1   ), p_arr( 2, 1   ), p_arr( 3, 1   ), 300, 'filled', 'linewidth', 3, 'markerfacecolor', [0 0.4470 0.7410], 'markeredgecolor', 'k'  );
    scatter3( anim.hAxes, p_arr( 1, end ), p_arr( 2, end ), p_arr( 3, end ), 300, 'filled', 'linewidth', 3, 'markerfacecolor', [0 0.4470 0.7410], 'markeredgecolor', 'k' );
    
    view([109.1242, 2.7291])

    for i = 1 : Ntmp
        robot.updateKinematics( q_arr( :, i ) );
        H = robot.getForwardKinematics( q_arr( :, i ) );
        anim.update( t_arr( i ) );
    end

    anim.close( );
end

%% Imitation Learning #1 Position 

% For Imitation Learning, one should define the number of basis function
N  = 100;

% Parameters of the 3 DMPs
alpha_z = 300.0;
alpha_s = 1.0;
beta_z  = 1/4 * alpha_z;
tau     = max( t_arr );

y0      = p_arr( :, 1   );
z0      = dp_arr( :, 1 );

% Get the number of total samples
[ ~, P_total ] = size( p_arr );

idx = 1 : 5 : P_total;
P = length( idx );

% Get the goal location  
g       = p_arr( :, idx( end ) );

cs      = CanonicalSystem( 'discrete', tau, alpha_s );
fs_x    = NonlinearForcingTerm( cs, N );
fs_y    = NonlinearForcingTerm( cs, N );
fs_z    = NonlinearForcingTerm( cs, N );

fs = cell( 1, 3 ); 
fs{1} = fs_x; fs{2} = fs_y; fs{3} = fs_z;

% Learning the N weights based on Locally Weighted Regression
w_arr = zeros( 3, N );

for k = 1 : 3
    for i = 1: N

        a_arr   = zeros( 1, P );
        f_arr   = zeros( 1, P );
        phi_arr = zeros( 1, P );

        for j = 1 : P
            a_arr( j ) = ( g( k ) - y0( k ) ) * cs.calc( t_arr( idx( j ) ) );

            ddy_des = data_raw.ddp_arr_filt( k, idx( j ) );
            dy_des  =  data_raw.dp_arr_filt( k, idx( j ) );
            y_des   =   data_raw.p_arr_filt( k, idx( j ) );

            % [ y_des, dy_des, ddy_des ] = min_jerk_traj( t_arr( idx( j ) ), y0( k ), g( k ), tau, 0 );

            f_arr( j ) = tau^2 * ddy_des + alpha_z * tau * dy_des + alpha_z * beta_z * ( y_des - g( k ) );
            phi_arr( j ) = fs{ k }.calc_ith( t_arr( idx( j ) ), i );
        end

        if( sum( a_arr .* phi_arr .* a_arr ) ~= 0 )
            w_arr( k, i ) = sum( a_arr .* phi_arr .* f_arr ) / sum( a_arr .* phi_arr .* a_arr );
        else
            w_arr( k, i ) = 0;
        end

    end

end

%% ---- (1C) Generating a Full trajectory with the Transformation System

trans_x = TransformationSystem( alpha_z, beta_z, tau, y0( 1 ), z0( 1 ) );
trans_y = TransformationSystem( alpha_z, beta_z, tau, y0( 2 ), z0( 2 ) );
trans_z = TransformationSystem( alpha_z, beta_z, tau, y0( 3 ), z0( 3 ) );
trans = cell( 1, 3 );
trans{ 1 } = trans_x; trans{ 2 } = trans_y; trans{ 3 } = trans_z;

% The time step of the simulation and its number of iteration
dt = t_arr( 2) - t_arr( 1 );
Nt = P_total;

% The total time and its time array
T  = dt * Nt;
t_arr2 = dt * (0:Nt);

% For plotting and saving the data
y_arr  = zeros( 3, Nt + 1 );
z_arr  = zeros( 3, Nt + 1 );
dz_arr = zeros( 3, Nt + 1 );

y_arr( :, 1 ) = y0;
z_arr( :, 1 ) = z0;

t = 0;
t0i = 0;

f_input_arr = zeros( 1, Nt );

for i = 1 : Nt
    
    for k = 1 : 3
        % Before conducting the movement
        % Maintaining that posture
        if t <= t0i
            y_arr( k, i + 1 ) = y0( k );
            z_arr( k, i + 1 ) = z0( k );

        % During the movement
        elseif t0i <= t

            % taking off the initial time offset
            t_tmp = t - t0i;

            % Calculating the input from the weights
            phi_act_arr = zeros( 1, N );
            for j = 1 : N
                phi_act_arr( j ) = fs{ k }.calc_ith( t_tmp, j );
            end

            if sum( phi_act_arr ) == 0
                f_input = 0;
            else
                f_input = sum( phi_act_arr .* w_arr( k, : ) )/sum( phi_act_arr ) * cs.calc( t_tmp ) * ( g( k ) - y0( k ) );
            end

            [ y, z, ~, dz ] = trans{ k }.step( g( k ), f_input, dt );
            y_arr(  k, i + 1 ) = y;
            z_arr(  k, i + 1 ) = z;
            dz_arr( k, i + 1 ) = dz;
        end

    end
    t = t + dt;
end

subplot( 3, 1, 1 )
hold on
plot( t_arr2, y_arr( 1, : ) )
plot( t_arr, data_raw.p_arr_filt( 1, : ) )


subplot( 3, 1, 2 )
hold on
plot( t_arr2, y_arr( 2, : ) )
plot( t_arr, data_raw.p_arr_filt( 2, : ) )

subplot( 3, 1, 3 )
hold on
plot( t_arr2, y_arr( 3, : ) )
plot( t_arr, data_raw.p_arr_filt( 3, : ) )

