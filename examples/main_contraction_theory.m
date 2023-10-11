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
D   = 3.0;
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
Nt = 8000;

% The total time and its time array
T     = dt * Nt;
t_arr = dt * (0:(Nt-1));

% For plotting and saving the data
y_arr  = zeros( 1, Nt );
dy_arr = zeros( 1, Nt );
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
        dy_arr( i + 1 ) = dy; 
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
plot( t_arr, dy_arr, 'linewidth', 3 )
plot(   t_P, dy_des_arr, 'linestyle', '--', 'linewidth', 5, 'color', 'k' )
set( gca, 'xlim', [0, T], 'fontsize', 30, 'xticklabel', {} )
ylabel( 'Vel. (-)' )

subplot( 3, 1, 3 )
hold on
plot( t_arr, dz_arr/D, 'linewidth', 3 )
plot(   t_P, ddy_des_arr, 'linestyle', '--', 'linewidth', 5, 'color', 'k' )
set( gca, 'xlim', [0, T], 'fontsize', 30, 'xtick', 0:0.5:T)
xlabel( 'Time (sec)' )
ylabel( 'Acc. (-)' )


%% [1D] Unrolling a plotting - Method 2
% From the learned weights, one can simply unroll a trajectory from
% state-space representation. 

% The time step of the simulation and its number of iteration
dt = 1e-3;
Nt = 8000;

% The total time and its time array
T     = dt * Nt;
t_arr = dt * (0:(Nt-1));

% First DMP
s_arr  = zeros( 1, Nt ); 
y_arr  = zeros( 2, Nt );
dy_arr = zeros( 2, Nt );

% Initial Position and Velocity 
y0  = [ 4.0, 2.0 ];
g   = [ 2.0, 3.0 ];
dy0 = [ 0.0, 0.0 ];

x_init = [ cs.calc( 0 ) ; y0'; dy0' ];
x_arr = zeros( 5, Nt );

Az = alpha_z*beta_z/tau^2*eye( 2 );
Bz =        alpha_z/tau*eye( 2 );

A1 = [   -alpha_s/tau, 0, 0, 0, 0;
        zeros( 2, 1 ), zeros(2, 2 ), eye( 2 );
        zeros( 2, 1 ),          -Az,    -Bz ];

x_curr = x_init;

t = 0;
t0i = 1.0;

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
        phi_sum1 = fs.calc_whole_at_t( t_tmp );
        
        f_input_x = 0;
        if phi_sum1 ~= 0
            f_input_x = fs.calc_whole_weighted_at_t( t_tmp, w_arr_LSS )/phi_sum1;
            f_input_x = f_input_x*( g( 1 )-y0( 1 ) )*x_curr( 1 );
        end

        phi_sum2 = fs.calc_whole_at_t( t_tmp );
        
        f_input_y = 0;
        if phi_sum2 ~= 0
            f_input_y = fs.calc_whole_weighted_at_t( t_tmp, w_arr_LSS )/phi_sum2;
            f_input_y = f_input_y*( g( 2 )-y0( 2 ) )*x_curr( 1 );
        end
        
        if t >= t0i + D
            f_input_x = 0;
            f_input_y = 0;
        end
        
        dx = A1 * x_curr + 1/tau^2*[ zeros(3,1); f_input_x;f_input_y ]+ [ zeros(3,1); Az*g' ];
        x_curr = x_curr + dx * dt;
        x_arr(:, i+1 ) = x_curr;
        
    end
    t = t + dt;
end

plot( t_arr, x_arr( 2:3, : ) )


%% [2A] Sequencing 2 Movements
% From the learned weights, one can simply unroll a trajectory from
% state-space representation. 

close all;

% The time step of the simulation and its number of iteration
dt = 1e-3;
Nt = 6000;

% The total time and its time array
T     = dt * Nt;
t_arr = dt * (0:(Nt-1));

% First DMP
s1_arr  = zeros( 1, Nt ); 
y1_arr  = zeros( 2, Nt );
dy1_arr = zeros( 2, Nt );

% Initial Position and Velocity of the First Movement
y0_1  = [  3.0; 0.0 ];
g1    = [  4.0; 2.0 ];
dy0_1 = [ 0.0; 0.0 ];

x1_init = [ cs.calc( 0 ); y0_1; dy0_1 ];
x1_arr = zeros( 5, Nt );

% Second DMP
s2_arr  = zeros( 1, Nt ); 
y2_arr  = zeros( 2, Nt );
dy2_arr = zeros( 2, Nt );

% Initial Position and Velocity of the Second Movement
y0_2  = [  0.0; 2.0 ];
g2    = [ 2.0; 0.0 ];
dy0_2 = [  0.0; 0.0 ];

% The DMPs
x2_init = [ cs.calc( 0 ) ; y0_2; dy0_2 ];
x2_arr = zeros( 5, Nt );

Az = alpha_z*beta_z/tau^2*eye( 2 );
Bz =        alpha_z/tau*eye( 2 );

A1 = [   -alpha_s/tau, 0, 0, 0, 0;
        zeros( 2, 1 ), zeros(2, 2 ), eye( 2 );
        zeros( 2, 1 ),          -Az,    -Bz ];

x1_curr  = x1_init;
x2_curr  = x2_init;

% Coupled
x12_curr = x2_init;
x12_arr = zeros( 5, Nt );

t = 0;
tmp_gain = 0;
t0i1 = 1.8;
t0i2 = 0.5;
% Forward Integration
for i = 0 : (Nt-1)

    % First DMP
    if t <= t0i1
        x1_arr( :, i+1 ) = x1_init;

    % During the movement
    elseif t0i1 <= t

        % taking off the initial time offset
        t_tmp = t - t0i1;

        % Calculating the input from the weights
        % First, check if whole activation value is 0
        phi_sum1 = fs.calc_whole_at_t( t_tmp );
        
        f_input_x = 0;
        if phi_sum1 ~= 0
            f_input_x = fs.calc_whole_weighted_at_t( t_tmp, w_arr_LSS )/phi_sum1;
            f_input_x = f_input_x*( g1( 1 )-y0_1( 1 ) )*x1_curr( 1 );
        end

        phi_sum2 = fs.calc_whole_at_t( t_tmp );
        
        f_input_y = 0;
        if phi_sum2 ~= 0
            f_input_y = fs.calc_whole_weighted_at_t( t_tmp, w_arr_LSS )/phi_sum2;
            f_input_y = f_input_y*( g1( 2 )-y0_1( 2 ) )*x1_curr( 1 );
        end

        % During the movement
        if t0i + D <= t        
            f_input_x = 0;
            f_input_y = 0;
        end

        dx = A1 * x1_curr + 1/tau^2*[ zeros(3,1); f_input_x;f_input_y ]+ [ zeros(3,1); Az*g1 ];
        x1_curr = x1_curr + dx * dt;
        x1_arr(:, i+1 ) = x1_curr;


    end
    
    % Second DMP
    if t <= t0i2
        x2_arr( :, i+1 ) = x2_init;

    % During the movement
    elseif t0i2 <= t

        % taking off the initial time offset
        t_tmp = t - t0i2;

        % Calculating the input from the weights
        % First, check if whole activation value is 0
        phi_sum1 = fs.calc_whole_at_t( t_tmp );
        
        f_input_x = 0;
        if phi_sum1 ~= 0
            f_input_x = fs.calc_whole_weighted_at_t( t_tmp, w_arr_LSS )/phi_sum1;
            f_input_x = f_input_x*( g2( 1 )-y0_2( 1 ) )*x2_curr( 1 );
        end

        phi_sum2 = fs.calc_whole_at_t( t_tmp );
        
        f_input_y = 0;
        if phi_sum2 ~= 0
            f_input_y = fs.calc_whole_weighted_at_t( t_tmp, w_arr_LSS )/phi_sum2;
            f_input_y = f_input_y*( g2( 2 )-y0_2( 2 ) )*x2_curr( 1 );
        end

        % During the movement
        if t0i + D <= t        
            f_input_x = 0;
            f_input_y = 0;
        end

        dx = A1 * x2_curr + 1/tau^2*[ zeros(3,1); f_input_x;f_input_y ]+ [ zeros(3,1); Az*g2 ];
        x2_curr = x2_curr + dx * dt;
        x2_arr( :, i+1 ) = x2_curr;

    end    


    % Coupled DMP
    if t <= t0i2
        x12_arr( :, i+1 ) = x2_init;

    % During the movement
    elseif t0i2 <= t

        % taking off the initial time offset
        t_tmp = t - t0i2;

        % Calculating the input from the weights
        % First, check if whole activation value is 0
        phi_sum1 = fs.calc_whole_at_t( t_tmp );
        
        f_input_x = 0;
        if phi_sum1 ~= 0
            f_input_x = fs.calc_whole_weighted_at_t( t_tmp, w_arr_LSS )/phi_sum1;
            f_input_x = f_input_x*( g2( 1 )-y0_2( 1 ) )*x12_curr( 1 );
        end

        phi_sum2 = fs.calc_whole_at_t( t_tmp );
        
        f_input_y = 0;
        if phi_sum2 ~= 0
            f_input_y = fs.calc_whole_weighted_at_t( t_tmp, w_arr_LSS )/phi_sum2;
            f_input_y = f_input_y*( g2( 2 )-y0_2( 2 ) )*x12_curr( 1 );
        end

        % Add on off 
        if t >= t0i2 + D*0.6;
            tmp_gain = tmp_gain + 0.004;
            if tmp_gain >= 1
               tmp_gain = 1; 
            end
            
        end

        dx = A1 * x12_curr + 1/tau^2*[ zeros(3,1); f_input_x;f_input_y ]+ [ zeros(3,1); Az*g2 ] + tmp_gain  * ( 2*eye( 5 ) * ( x1_curr - x12_curr ) + [ zeros(3,1); Az*g1 ] - [ zeros(3,1); Az*g2 ] );        
        
        x12_curr = x12_curr + dx * dt;
        x12_arr(:, i+1 ) = x12_curr;

    end    
    
    
    t = t + dt;
end

% Draw an Animation 
f = figure( ); a = axes( 'parent', f );
hold on
set( a, 'xlim', [-1, 3], 'ylim', [-1,3 ], 'fontsize', 30 );
axis equal
% First DMP
xlabel( 'X (m)', 'fontsize', 35 );
ylabel( 'Y (m)', 'fontsize', 35 );
scatter( a, y0_1( 1 ), y0_1( 2 ), 300, 'o', 'filled', 'markerfacecolor', [0, 0.4470, 0.7410], 'markeredgecolor', 'black', 'markerfacealpha', 0.8 )
scatter( a,  g1( 1 ),   g1( 2 ), 300,'square', 'filled', 'markerfacecolor', [0, 0.4470, 0.7410], 'markeredgecolor', 'black', 'markerfacealpha', 0.8 )
plot( a, x1_arr( 2, : ), x1_arr( 3, : ), 'linewidth', 3, 'color', 'black' )

pDMP1 = scatter( a, x1_arr( 2, 1 ), x1_arr( 3, 1 ), 500, 'o', 'filled', 'markerfacecolor', [0, 0.4470, 0.7410], 'markeredgecolor', 'black', 'markerfacealpha', 0.8 );

% Second DMP
scatter( a, y0_2( 1 ), y0_2( 2 ), 300, 'o', 'filled', 'markerfacecolor', [0.8500, 0.3250, 0.0980], 'markeredgecolor', 'black', 'markerfacealpha', 0.8 )
scatter( a,  g2( 1 ),   g2( 2 ), 300,'square', 'filled', 'markerfacecolor', [0.8500, 0.3250, 0.0980], 'markeredgecolor', 'black', 'markerfacealpha', 0.8 )
plot( a, x2_arr( 2, : ), x2_arr( 3, : ), 'linewidth', 3, 'color', 'black' )

pDMP2 = scatter( a, x2_arr( 2, 1 ), x2_arr( 3, 1 ), 500, 'o', 'filled', 'markerfacecolor', [0.8500, 0.3250, 0.0980], 'markeredgecolor', 'black', 'markerfacealpha', 0.8 );

pDMP12  = scatter( a, x12_arr( 2, 1 ), x12_arr( 3, 1 ), 500, 'o', 'filled', 'markerfacecolor', [0.4940, 0.1840, 0.5560]	, 'markeredgecolor', 'black', 'markerfacealpha', 0.8 );

tmp_step = round( 1/dt / 30 );

v = VideoWriter( 'video.mp4','MPEG-4' );
v.FrameRate = 30;
open( v );

for i = 1 : tmp_step : Nt
    set( pDMP1, 'XData', x1_arr( 2, i ), 'YData', x1_arr( 3, i ) );
    set( pDMP2, 'XData', x2_arr( 2, i ), 'YData', x2_arr( 3, i ) );
    
    set( pDMP12, 'XData', x12_arr( 2, i ), 'YData', x12_arr( 3, i ) );
    drawnow 
    i
    
    tmp_frame = getframe( f );
    writeVideo( v,tmp_frame );

end
close( v );

%% [2B] Sequencing 8 Movements

close all;
% The number of movements 
Nmov = 8;

% The time step of the simulation and its number of iteration
dt = 1e-3;
Nt = 3500 * Nmov;

% The total time and its time array
T     = dt * Nt;
t_arr = dt * (0:(Nt-1));

% The Nmov DMPs
s_arr  = zeros( 1, Nt, Nmov );
y_arr  = zeros( 2, Nt, Nmov );
dy_arr = zeros( 2, Nt, Nmov );
x_arr  = zeros( 5, Nt, Nmov );

x_curr = zeros( 5, Nmov );

x_coupled = zeros( 5, 1 );

% Initial Position and Velocity of the First Movement
y0 = zeros( 2, Nmov );
g  = zeros( 2, Nmov );
g( :, 1 ) = [ 2, 4 ];


for i = 2 : Nmov
   y0( :, i ) = g( :, i-1 ) + [4;0];
    g( :, i ) = y0( :, i ) + [ g( 1, 1 ); 2*(mod( i, 2 )-0.5 )*g( 2, 1 ) ] ;
end

for i = 1 : Nmov
   x_curr( :, i ) = [ cs.calc( 0 ); y0( :, i ); zeros( 2, 1 ) ];
end

x_coupled = x_curr( :, 1 );
x_coupled_arr = zeros( 5, Nmov );
x_coupled_arr( :, 1 ) = x_coupled;

x_arr( :, 1, : ) = x_curr;

f = figure( ); a = axes( 'parent', f );
hold on 
axis equal
set( a, 'xlim', [-5, max( y0( 1, : ) )+5 ], 'ylim', [-1, max( y0( 2, : ) )+1 ] )

for i = 1 : Nmov
   scatter( a, y0( 1, i ), y0( 2, i ), 200,      'o', 'filled', 'markerfacecolor', [0.8500, 0.3250, 0.0980], 'markeredgecolor', 'black', 'markerfacealpha', 0.8 );
   scatter( a,  g( 1, i ),  g( 2, i ), 200, 'square', 'filled', 'markerfacecolor', [0.8500, 0.3250, 0.0980], 'markeredgecolor', 'black', 'markerfacealpha', 0.8 );
end

Az = alpha_z*beta_z/tau^2*eye( 2 );
Bz =        alpha_z/tau*eye( 2 );

A1 = [   -alpha_s/tau, 0, 0, 0, 0;
        zeros( 2, 1 ), zeros(2, 2 ), eye( 2 );
        zeros( 2, 1 ),          -Az,    -Bz ];
    
t = 0;

t0i   =  0.5;
toff  =  -1.0;
toff2 =  0;

acts = zeros( 1, Nmov-1 );

for i = 0 : (Nt-1)

    for j = 1 : Nmov
        
       if t >= (t0i+(j-1)*(D+toff))
           
            % taking off the initial time offset
            t_tmp = t - (t0i+(j-1)*(D+toff));

            % Calculating the input from the weights
            % First, check if whole activation value is 0
            phi_sum1 = fs.calc_whole_at_t( t_tmp );

            f_input_x = 0;
            if phi_sum1 ~= 0
                f_input_x = fs.calc_whole_weighted_at_t( t_tmp, w_arr_LSS )/phi_sum1;
                f_input_x = f_input_x*( g( 1, j )-y0( 1, j ) )*x_curr( 1, j );
            end

            phi_sum2 = fs.calc_whole_at_t( t_tmp );

            f_input_y = 0;
            if phi_sum2 ~= 0
                f_input_y = fs.calc_whole_weighted_at_t( t_tmp, w_arr_LSS )/phi_sum2;
                f_input_y = f_input_y*( g( 2, j )-y0( 2, j ) )*x_curr( 1, j );
            end
            
            if t >=  (t0i+(j-1)*(D+toff)+D)
                f_input_x = 0;
                f_input_y = 0;
            end
            
            dx = A1 * x_curr( :, j ) + 1/tau^2*[ zeros(3,1); f_input_x;f_input_y ]+ [ zeros(3,1); Az*g( :, j) ];

            x_curr( :, j ) = x_curr( :, j ) + dx * dt;
            x_arr( :, i+1, j ) = x_curr( :, j );           
           
       else    
           x_arr( :, i+1, j ) = [ cs.calc( 0 ); y0( :, j ); zeros( 2, 1 ) ];
       end
       
    end
    
    % For the coupled movements
    if t >= t0i
        for j = 1 : Nmov

           if t >= (t0i+(j-1)*(D+toff))

                % taking off the initial time offset
                t_tmp = t - (t0i+(j-1)*(D+toff));

                % Calculating the input from the weights
                % First, check if whole activation value is 0
                phi_sum1 = fs.calc_whole_at_t( t_tmp );

                f_input_x = 0;
                if phi_sum1 ~= 0
                    f_input_x = fs.calc_whole_weighted_at_t( t_tmp, w_arr_LSS )/phi_sum1;
                    f_input_x = f_input_x*( g( 1, j )-y0( 1, j ) )*x_curr( 1, j );
                end

                phi_sum2 = fs.calc_whole_at_t( t_tmp );

                f_input_y = 0;
                if phi_sum2 ~= 0
                    f_input_y = fs.calc_whole_weighted_at_t( t_tmp, w_arr_LSS )/phi_sum2;
                    f_input_y = f_input_y*( g( 2, j )-y0( 2, j ) )*x_curr( 1, j );
                end

                if t >= (t0i+(j-1)*(D+toff)+D)
                    f_input_x = 0;
                    f_input_y = 0;
                end

                dx_coupled = A1 * x_coupled + 1/tau^2*[ zeros(3,1); f_input_x;f_input_y ]+ [ zeros(3,1); Az*g( :, j) ];

                % Need to calculate the the coupling term 
                if  t >= (t0i+(j-1)*D + 0.5*D+toff2 ) && j ~= 8
                    acts( j ) = acts( j ) + 0.001;

                    if acts( j ) >= 1.0
                       acts( j ) = 1.0; 
                    end
                    dx_coupled = dx_coupled + acts( j ) * ( 2 * eye( 5 ) * ( x_arr( :, i, j + 1 ) -x_arr( :, i, j ) ) + [ zeros(3,1); Az*g( :, j+1) ] - [ zeros(3,1); Az*g( :, j) ] );
                end

                x_coupled = x_coupled + dx_coupled * dt;
                x_coupled_arr( :, i+1 ) = x_coupled;           
           end
        end
    else
        x_coupled_arr( :, i+1 ) = [ cs.calc( 0 ); y0( :, 1 ); zeros( 2, 1 ) ];
    end    
    
    
    t = t + dt;
end

pTrack = cell( 1, 8 );

for i = 1 : Nmov
    plot( a, x_arr( 2, :, i ), x_arr( 3, :, i ), 'linewidth', 3, 'color', 'black' )
    pTrack{ i } = scatter( a, x_arr( 2, 1, i ), x_arr( 3, 1, i ), 400, 'filled', ...
                        'markerfacealpha', 0.7, 'markerfacecolor', [0.4940, 0.1840, 0.5560], 'markeredgecolor', 'black' );
end

pOrig = scatter( a, x_coupled_arr( 2, 1 ), x_coupled_arr( 3, 1 ), 800, 'filled', ...
                    'markerfacealpha', 0.7, 'markerfacecolor', [0.4940, 0.1840, 0.5560], 'markeredgecolor', 'black' );


v = VideoWriter( 'video_whole.mp4','MPEG-4' );
v.FrameRate = 30;
open( v );
tmp_step = 33;
for i = 1 : tmp_step : Nt
    
    for j = 1 : 8
        set( pTrack{ j }, 'XData', x_arr( 2, i, j ), 'YData', x_arr( 3, i, j ) );
    end
    
    set( pOrig, 'XData', x_coupled_arr( 2, i ), 'YData', x_coupled_arr( 3, i ) ); 
    drawnow 
    
    tmp_frame = getframe( f );
    writeVideo( v,tmp_frame );
    i
end
close( v );
