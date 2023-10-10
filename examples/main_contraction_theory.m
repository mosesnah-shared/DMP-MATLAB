%% [Example Script] Contraction Theory and DMP
% [Author] Moses Chong-ook Nah
% [Date]   2023.10.08

%% [Initialization] 
clear; close all; clc;
fig_config( 'fontSize', 20, 'markerSize', 10 )

%% [1A] [Training A Single 2D Movement] 
% There is two transformation system
% In 2D Space, we aim to imitate the Minimum jerk Trajectory
% The initial (p0i), final posture (p0f), duration (D), starting time (t0i) of the trajectory
p0i = 0.0;
p0f = 1.0;
D   = 1.0;
t0i = 0.0;

% For Imitation Learning, one should define the number of basis function
N  = 50;

% Parameters of the 3 DMPs
alpha_z = 10.0;
alpha_s = 1.0;
beta_z  = 1/4 * alpha_z;
tau     = D;
y0      = p0i;
g       = p0f;
z0      = 0;

% A single canonical system
cs        = CanonicalSystem( 'discrete', tau, alpha_s );
trans_sys = TransformationSystem( alpha_z, beta_z, tau, y0, z0 );
fs        = NonlinearForcingTerm( cs, N );

%% [1B] [Learning Weights via Locally Weighted Regression or Least-Square]

% Assume that we have P sample points for the minimum jerk trajectory.
P = 100;

% Equal Sampling along time with duration D
t_P = linspace( 0.0, D, P );

% Desired Trajectories
  y_des_arr = zeros( 1, P );
 dy_des_arr = zeros( 1, P );
ddy_des_arr = zeros( 1, P );

for i = 1 : P
    [ y_des, dy_des, ddy_des ] = min_jerk_traj( t_P( i ), p0i, p0f, D, 0 );
      y_des_arr( i ) =   y_des;
     dy_des_arr( i ) =  dy_des;
    ddy_des_arr( i ) = ddy_des;
end

% Calculating a_arr and f_arr that can be calculated
f_arr = trans_sys.get_desired( y_des_arr, dy_des_arr, ddy_des_arr, g );

w_arr_LSS = zeros( 1, N );

phi_mat = zeros( P, N );
for i = 1 : P
    phi_sum = fs.calc_whole_at_t( t_P( i ) );
    for j = 1 : N
        phi_mat( i, j ) = fs.calc_ith( t_P( i ), j ) / phi_sum * ( g - y0 ) * cs.calc( t_P( i ) );
    end
end

% Get w_arr with Least square solution
tmp = ( phi_mat' * phi_mat )^(-1) * phi_mat' * f_arr';
w_arr_LSS( : ) = tmp';

%% [1C] Unrolling a plotting - Method 1
trans_sys.reset( )

% The time step of the simulation and its number of iteration
dt = 1e-3;
Nt = 2000;

% The total time and its time array
T     = dt * Nt;
t_arr = dt * (0:(Nt-1));

% For plotting and saving the data
y_arr  = zeros( 1, Nt );
z_arr  = zeros( 1, Nt );
dz_arr = zeros( 1, Nt );

y_arr( 1 ) = y0;
z_arr( 1 ) = z0;

t = 0;
f_input_arr = zeros( 1, Nt );

% One can also check spatial temporance from the weight
% for tau = 0.5:0.1:2.0
% for   g = 0.5:0.1:2.0

for i = 0 : (Nt-1)
    
    % Before conducting the movement
    % Maintaining that posture
    if t <= t0i
        y_arr( i + 1 ) = y0;
        z_arr( i + 1 ) = z0;
        
    % During the movement
    elseif t0i <= t
        
        if t<= t0i + D

            % taking off the initial time offset
            t_tmp = t - t0i;

            % Calculating the input from the weights
            % First, check if whole activation value is 0
            phi_sum = fs.calc_whole_at_t( t_tmp );
            
            f_input = 0;

            if phi_sum ~= 0
                f_input = fs.calc_whole_weighted_at_t( t_tmp, w_arr_LSS )/phi_sum;
                f_input = f_input*(g-y0)*cs.calc( t_tmp );
            end
            
        else
            f_input = 0; 
        end

        f_input_arr( i ) = f_input; 
        [ y, z, dy, dz ] = trans_sys.step( g, f_input, dt );
        y_arr(  i + 1 ) = y;
        z_arr(  i + 1 ) = z;
        dz_arr( i + 1 ) = dz;        
        dz_arr( i + 1 ) = dz;
    end

    t = t + dt;
end

subplot( 3, 1, 1 )
hold on
plot( t_arr, y_arr, 'linewidth', 3 )
plot( t_P, y_des_arr, 'linestyle', '--', 'linewidth', 5, 'color', 'k' )
set( gca, 'xlim', [0, T], 'fontsize', 30, 'xticklabel', {} )
ylabel( 'Pos. (-)' )

subplot( 3, 1, 2 )
hold on
plot( t_arr, z_arr, 'linewidth', 3 )
plot(   t_P, dy_des_arr, 'linestyle', '--', 'linewidth', 5, 'color', 'k' )
set( gca, 'xlim', [0, T], 'fontsize', 30, 'xticklabel', {} )
ylabel( 'Vel. (-)' )

subplot( 3, 1, 3 )
hold on
plot( t_arr, dz_arr, 'linewidth', 3 )
plot(   t_P, ddy_des_arr, 'linestyle', '--', 'linewidth', 5, 'color', 'k' )
set( gca, 'xlim', [0, T], 'fontsize', 30, 'xtick', 0:0.5:T)
xlabel( 'Time (sec)' )
ylabel( 'Acc. (-)' )


%% [1C] Unrolling a plotting - Method 2
% From the learned weights, one can simply unroll a trajectory from
% state-space representation. 

% The time step of the simulation and its number of iteration
dt = 1e-3;
Nt = 2000;

% The total time and its time array
T     = dt * Nt;
t_arr = dt * (0:(Nt-1));

% First DMP
s_arr  = zeros( 1, Nt ); 
y_arr  = zeros( 2, Nt );
dy_arr = zeros( 2, Nt );

% Initial Position and Velocity 
y0  = [ 1.0, 2.0 ];
g   = [ 4.0, 1.0 ];
dy0 = [ 0.0, 0.0 ];

cs   = CanonicalSystem( 'discrete', tau, alpha_s );
fs_x = NonlinearForcingTerm( cs, N );
fs_y = NonlinearForcingTerm( cs, N );

x_init = [ cs.calc( 0 ) ; y0'; dy0' ];
x_arr = zeros( 5, Nt );

Az = alpha_z*beta_z/tau^2*eye( 2 );
Bz =        alpha_z/tau*eye( 2 );

A1 = [   -alpha_s/tau, 0, 0, 0, 0;
        zeros( 2, 1 ), zeros(2, 2 ), eye( 2 );
        zeros( 2, 1 ),          -Az,    -Bz ];

x_curr = x_init;

t = 0;

% Forward Integration
for i = 0 : (Nt-1)

    % First DMP
    if t <= t0i
        x_arr( :, i+1 ) = x_init;

    % During the movement
    elseif t0i <= t

        % taking off the initial time offset
        t_tmp = t - t0i;

        % Calculating the input from the weights
        % First, check if whole activation value is 0
        phi_sum1 = fs_x.calc_whole_at_t( t_tmp );
        
        f_input_x = 0;
        if phi_sum1 ~= 0
            f_input_x = fs_x.calc_whole_weighted_at_t( t_tmp, w_arr_LSS )/phi_sum1;
            f_input_x = f_input_x*( g( 1 )-y0( 1 ) )*x_curr( 1 );
        end

        phi_sum2 = fs_y.calc_whole_at_t( t_tmp );
        
        f_input_y = 0;
        if phi_sum2 ~= 0
            f_input_y = fs_y.calc_whole_weighted_at_t( t_tmp, w_arr_LSS )/phi_sum2;
            f_input_y = f_input_y*( g( 2 )-y0( 2 ) )*x_curr( 1 );
        end

        dx = A1 * x_curr + 1/tau*[ zeros(3,1); f_input_x;f_input_y ]+ [ zeros(3,1); Az*g' ];
        x_curr = x_curr + dx * dt;
        x_arr(:, i+1 ) = x_curr;

        f_input_x
    end
    t = t + dt;
end

plot( t_arr, x_arr( 2:3, : ) )


%% [2A] Sequencing
% From the learned weights, one can simply unroll a trajectory from
% state-space representation. 

% The time step of the simulation and its number of iteration
dt = 1e-3;
Nt = 8000;

% The total time and its time array
T     = dt * Nt;
t_arr = dt * (0:(Nt-1));

% First DMP
s1_arr  = zeros( 1, Nt ); 
y1_arr  = zeros( 2, Nt );
dy1_arr = zeros( 2, Nt );

% Initial Position and Velocity of the First Movement
y0_1  = [ 1.0; 2.0 ];
g1   = [ 3.0; 1.0 ];
dy0_1 = [ 0.0; 0.0 ];

cs1 = CanonicalSystem( 'discrete', tau, alpha_s );
fs1_x = NonlinearForcingTerm( cs1, N );
fs1_y = NonlinearForcingTerm( cs1, N );

x1_init = [ cs1.calc( 0 ) ; y0_1; dy0_1 ];
x1_arr = zeros( 5, Nt );

% First DMP
s2_arr  = zeros( 1, Nt ); 
y2_arr  = zeros( 2, Nt );
dy2_arr = zeros( 2, Nt );

% Initial Position and Velocity of the First Movement
y0_2  = [ -0.4; 2.0 ];
g2   = [ 0.6; 1.5 ];
dy0_2 = [ 0.0; 0.0 ];

cs2 = CanonicalSystem( 'discrete', tau, alpha_s );
fs2_x = NonlinearForcingTerm( cs2, N );
fs2_y = NonlinearForcingTerm( cs2, N );

x2_init = [ cs2.calc( 0 ) ; y0_2; dy0_2 ];
x2_arr = zeros( 5, Nt );

Az = alpha_z*beta_z/tau^2*eye( 2 );
Bz =        alpha_z/tau*eye( 2 );

A1 = [   -alpha_s/tau, 0, 0, 0, 0;
        zeros( 2, 1 ), zeros(2, 2 ), eye( 2 );
        zeros( 2, 1 ),          -Az,    -Bz ];

x1_curr = x1_init;
x2_curr = x2_init;

t = 0;

% Forward Integration
for i = 0 : (Nt-1)

    % First DMP
    if t <= t0i +0.3
        x1_arr( :, i+1 ) = x1_init;

    % During the movement
    elseif t0i +0.3 <= t

        % taking off the initial time offset
        t_tmp = t - t0i;

        % Calculating the input from the weights
        % First, check if whole activation value is 0
        phi_sum1 = fs1_x.calc_whole_at_t( t_tmp );
        
        f_input_x = 0;
        if phi_sum1 ~= 0
            f_input_x = fs1_x.calc_whole_weighted_at_t( t_tmp, w_arr_LSS )/phi_sum1;
            f_input_x = f_input_x*( g1( 1 )-y0_1( 1 ) )*x1_curr( 1 );
        end

        phi_sum2 = fs1_y.calc_whole_at_t( t_tmp );
        
        f_input_y = 0;
        if phi_sum2 ~= 0
            f_input_y = fs1_y.calc_whole_weighted_at_t( t_tmp, w_arr_LSS )/phi_sum2;
            f_input_y = f_input_y*( g1( 2 )-y0_1( 2 ) )*x1_curr( 1 );
        end

        % During the movement
        if t0i + D <= t        
            f_input_x = 0;
            f_input_y = 0;
        end

        dx = A1 * x1_curr + 1/tau*[ zeros(3,1); f_input_x;f_input_y ]+ [ zeros(3,1); Az*g1 ];
        x1_curr = x1_curr + dx * dt;
        x1_arr(:, i+1 ) = x1_curr;


    end

    % First DMP
    if t <= t0i
        x2_arr( :, i+1 ) = x2_init;

    % During the movement
    elseif t0i <= t

        % taking off the initial time offset
        t_tmp = t - t0i;

        % Calculating the input from the weights
        % First, check if whole activation value is 0
        phi_sum1 = fs2_x.calc_whole_at_t( t_tmp );
        
        f_input_x = 0;
        if phi_sum1 ~= 0
            f_input_x = fs2_x.calc_whole_weighted_at_t( t_tmp, w_arr_LSS )/phi_sum1;
            f_input_x = f_input_x*( g2( 1 )-y0_2( 1 ) )*x2_curr( 1 );
        end

        phi_sum2 = fs2_y.calc_whole_at_t( t_tmp );
        
        f_input_y = 0;
        if phi_sum2 ~= 0
            f_input_y = fs2_y.calc_whole_weighted_at_t( t_tmp, w_arr_LSS )/phi_sum2;
            f_input_y = f_input_y*( g2( 2 )-y0_2( 2 ) )*x2_curr( 1 );
        end

        dx = A1 * x2_curr + 1/tau*[ zeros(3,1); f_input_x;f_input_y ]+ [ zeros(3,1); Az*g2 ] +  4 * eye( 5 ) * ( x1_curr - x2_curr ) + + [ zeros(3,1); Az*g1 ] - + [ zeros(3,1); Az*g2 ] ;
        x2_curr = x2_curr + dx * dt;
        x2_arr(:, i+1 ) = x2_curr;

    end

    t = t + dt;
end

plot( x1_arr( 2, : ), x1_arr( 3, : ) )
hold on
scatter( g1( 1 ), g1( 2 ), 400, "yellow", "filled" )
plot( x2_arr( 2, : ), x2_arr( 3, : ) )
scatter( g2( 1 ), g2( 2 ), 400, "green" , "filled")
axis equal