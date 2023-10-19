%% [Example Script] Contraction Theory and DMP
% [Author] Moses Chong-ook Nah
% [Date]   2023.10.08

%% [Initialization] 
clear; close all; clc;
fig_config( 'fontSize', 20, 'markerSize', 10 )

%% [1A] [Calling a Weights]

load( 'learned_parameters/min_jerk.mat' );
data1 = data;
data4 = data;

% Load the dataset 
load( 'learned_parameters/cosine.mat'   );
data2 = data;

load( 'learned_parameters/radial.mat'   );
data3 = data;

data_whole = { data1, data2, data3, data4 };

data_whole{ 1 }.p0i = [ 3, 3, 3 ];
data_whole{ 1 }.p0f = data_whole{ 1 }.p0i + data_whole{ 1 }.p0f;

data_whole{ 2 }.p0i = data_whole{ 1 }.p0f + [ 1, 1, 0 ];
data_whole{ 2 }.p0f = data_whole{ 2 }.p0i + data_whole{ 2 }.p0f;

data_whole{ 3 }.p0i;
data_whole{ 3 }.p0f = data_whole{ 3 }.p0i + data_whole{ 3 }.p0f;

data_whole{ 4 }.p0i = data_whole{ 3 }.p0f + [ 1,2, 1];
data_whole{ 4 }.p0f = data_whole{ 4 }.p0i + data_whole{ 4 }.p0f;


%% [1B] Sequencing the Movements

close all;

% The time step of the simulation and its number of iteration
dt = 1e-3;
Nt = 20000;

% The total time and its time array
T     = dt * Nt;
t_arr = dt * (0:(Nt-1));

% There are Four movements, in the following order
% (1) Min-jerk-trajectory
% (2) M-shape Looking Trajectory
% (3) Swirl-like Trajectory
% (4) Min-jerk-trajectory

% The x_arr of the 4 movements
% The state-vector is 1 + 3 + 3 = 7
x_arr  = zeros( 4, 7, Nt );
dx_arr = zeros( 4, 7 );
x0_arr = zeros( 4, 7 );

% The State-array 
Az_arr = zeros( 4, 3, 3 );
Bz_arr = zeros( 4, 3, 3 );
A_arr  = zeros( 4, 7, 7 );

for i = 1 : 4
    az  = data_whole{ i }.alpha_z;
    bz  = data_whole{ i }.beta_z;
    tau = data_whole{ i }.tau;
    
    Az = az*bz/tau^2*eye( 3 );
    Bz =    az/tau^1*eye( 3 );

    Az_arr( i, :, : ) = Az;
    Bz_arr( i, :, : ) = Bz;

    alpha_s = data_whole{ i }.alpha_s;
    
    A_arr( i, :, : ) = [   -alpha_s/tau,   zeros( 1, 3 ),   zeros( 1, 3 );
                          zeros( 3, 1 ),   zeros( 3, 3 ),        eye( 3 );
                          zeros( 3, 1 ),             -Az,           -Bz ];
    
end

% Setting the Initial Condition
for i = 1 : 4
    x_arr( i, 1  , 1 ) = data_whole{ i }.cs.calc( 0 );
    x_arr( i, 2:4, 1 ) = data_whole{ i }.p0i ;
    x_arr( i, 5:7, 1 ) = double( subs( data_whole{ i }.dp_sym, 0 ) ) ;

    x0_arr( i, : ) = x_arr( i, :, 1 );
end

% The initial time start 
t0i_arr = [ 0, 3.0, 9.0, 15.0 ];

x_curr = x0_arr;
t = 0;
% Forward Integration
for i = 0 : (Nt-1)

    % Integrating over the movement
    for j = 1 : 4
        
        if t <= t0i_arr( j )
            x_arr( j, :, i + 1 ) = x0_arr( j, : );
            x_curr( j, : ) = x0_arr( j, : );

        else
            % taking off the initial time offset
            t_tmp = t - t0i_arr( j );

            % Getting the 3D nonlinear forcing term 
            f_arr = zeros( 3, 1 );

            y0    = data_whole{ j }.p0i;                
            g     = data_whole{ j }.p0f;
            w_arr = data_whole{ j }.weight;

            for k = 1 : 3
                % Calculating the input from the weights
                % First, check if whole activation value is 0
                phi_sum = data_whole{ j }.fs.calc_whole_at_t( t_tmp );

                if phi_sum ~= 0
                    f_arr( k ) = data_whole{ j }.fs.calc_whole_weighted_at_t( t_tmp, w_arr( k, : ) )/phi_sum;
                    f_arr( k ) = f_arr( k )*( g( k )-y0( k ) )*x_curr( j, 1 );
                end

            end

            % After the movement
            if t0i_arr( j ) + data_whole{ j }.tau <= t        
                f_arr = zeros( 3, 1 );
            end         

            % Update the State Vector
            A  = squeeze(  A_arr( j, :, : ) );
            Az = squeeze( Az_arr( j, :, : ) );
    
            dx = A * x_curr( j, : )' + 1/data_whole{ j }.tau^2*[ zeros(4,1); f_arr ]+ [ zeros(4,1); Az*data_whole{ j }.p0f' ];
            x_curr( j, : ) = x_curr( j, : ) + dx' * dt;
            x_arr( j, :, i+1 ) = x_curr( j, : );            

        end
      
    end
    t = t+dt;
end

%% Plotting the Figures 
x1 = squeeze( x_arr( 1, 2, : ) );
y1 = squeeze( x_arr( 1, 3, : ) );
z1 = squeeze( x_arr( 1, 4, : ) );

x2 = squeeze( x_arr( 2, 2, : ) );
y2 = squeeze( x_arr( 2, 3, : ) );
z2 = squeeze( x_arr( 2, 4, : ) );

x3 = squeeze( x_arr( 3, 2, : ) );
y3 = squeeze( x_arr( 3, 3, : ) );
z3 = squeeze( x_arr( 3, 4, : ) );

x4 = squeeze( x_arr( 4, 2, : ) );
y4 = squeeze( x_arr( 4, 3, : ) );
z4 = squeeze( x_arr( 4, 4, : ) );

f = figure( ); a = axes( 'parent', f );
hold on;

plot3( x1, y1, z1, 'color', 'black', 'linewidth', 3 )
plot3( x2, y2, z2, 'color', 'black', 'linewidth', 3 )
plot3( x3, y3, z3, 'color', 'black', 'linewidth', 3 )
plot3( x4, y4, z4, 'color', 'black', 'linewidth', 3 )

