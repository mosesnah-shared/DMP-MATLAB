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

data_whole{ 1 }.p0i = [ 1, 1, 1 ];
data_whole{ 1 }.p0f = data_whole{ 1 }.p0i + data_whole{ 1 }.p0f;

data_whole{ 2 }.p0i = [ 5, 2, 2 ];
data_whole{ 2 }.p0f = data_whole{ 2 }.p0i + data_whole{ 2 }.p0f;

data_whole{ 3 }.p0i = [ 2, 4, 2 ];
data_whole{ 3 }.p0f = data_whole{ 3 }.p0i + data_whole{ 3 }.p0f;

data_whole{ 4 }.p0i = [ 1, 4, 4 ];
data_whole{ 4 }.p0f = data_whole{ 4 }.p0i + data_whole{ 4 }.p0f;


%% [1B] Sequencing the Movements

close all;

% The time step of the simulation and its number of iteration
dt = 1e-3;
Nt = 14000;

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
t0i_arr = [ 0, 3.0, 7.0, 12.0 ];

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


%%
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
Nt = 18500;

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
acts = zeros( 1, Nmov-1 );
k = 1;
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
        % Adding the force term 
        if k <= 8
            if t >= (t0i+(k-1)*(D+toff)) && t <= (t0i+k*(D+toff)) 
    
                % taking off the initial time offset
                t_tmp = t - (t0i+(k-1)*(D+toff));
    
                % Calculating the input from the weights
                % First, check if whole activation value is 0
                phi_sum1 = fs.calc_whole_at_t( t_tmp );
    
                f_input_x = 0;
                if phi_sum1 ~= 0
                    f_input_x = fs.calc_whole_weighted_at_t( t_tmp, w_arr_LSS )/phi_sum1;
                    f_input_x = f_input_x*( g( 1, k )-y0( 1, k ) )*x_curr( 1, k );
                end
    
                phi_sum2 = fs.calc_whole_at_t( t_tmp );
    
                f_input_y = 0;
                if phi_sum2 ~= 0
                    f_input_y = fs.calc_whole_weighted_at_t( t_tmp, w_arr_LSS )/phi_sum2;
                    f_input_y = f_input_y*( g( 2, k )-y0( 2, k ) )*x_curr( 1, k );
                end
                
                if t >=  (t0i+(k-1)*(D+toff)+D)
                    f_input_x = 0;
                    f_input_y = 0;
                end
                
                dx_coupled = A1 * x_coupled + 1/tau^2*[ zeros(3,1); f_input_x;f_input_y ]+ [ zeros(3,1); Az*g( :, k) ];
               
                    
                dx_coupled = dx_coupled + ( 2*eye( 5 ) * ( x_curr( :, k ) - x_coupled ) ); 
    
                x_coupled = x_coupled + dx_coupled * dt;
                x_coupled_arr( :, i+1 ) = x_coupled;    
                
            else
                k = k+1;
            end
        else

                % taking off the initial time offset
                t_tmp = t - (t0i+(8-1)*(D+toff));
    
                % Calculating the input from the weights
                % First, check if whole activation value is 0
                phi_sum1 = fs.calc_whole_at_t( t_tmp );
    
                f_input_x = 0;
                if phi_sum1 ~= 0
                    f_input_x = fs.calc_whole_weighted_at_t( t_tmp, w_arr_LSS )/phi_sum1;
                    f_input_x = f_input_x*( g( 1, 8 )-y0( 1, 8 ) )*x_curr( 1, 8 );
                end
    
                phi_sum2 = fs.calc_whole_at_t( t_tmp );
    
                f_input_y = 0;
                if phi_sum2 ~= 0
                    f_input_y = fs.calc_whole_weighted_at_t( t_tmp, w_arr_LSS )/phi_sum2;
                    f_input_y = f_input_y*( g( 2, 8 )-y0( 2, 8 ) )*x_curr( 1, 8 );
                end
                
                if t >=  (t0i+(8-1)*(D+toff)+D)
                    f_input_x = 0;
                    f_input_y = 0;
                end
                
                dx_coupled = A1 * x_coupled + 1/tau^2*[ zeros(3,1); f_input_x;f_input_y ]+ [ zeros(3,1); Az*g( :, 8) ];
               
                    
                dx_coupled = dx_coupled + ( 2*eye( 5 ) * ( x_curr( :, 8 ) - x_coupled ) ); 
    
                x_coupled = x_coupled + dx_coupled * dt;
                x_coupled_arr( :, i+1 ) = x_coupled;                

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
