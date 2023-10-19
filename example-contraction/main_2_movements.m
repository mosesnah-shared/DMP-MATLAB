%% [Example Script] Contraction Theory and DMP
% [Author] Moses Chong-ook Nah
% [Date]   2023.10.08

%% [Initialization] 
clear; close all; clc;
fig_config( 'fontSize', 20, 'markerSize', 10 )

%% [1A] [Calling a Weights]

% Load the dataset 
load( 'learned_parameters/cosine.mat'   );
data1 = data;

load( 'learned_parameters/min_jerk.mat' );
data2 = data;

load( 'learned_parameters/radial.mat'   );
data3 = data;

%% [1B] Sequencing 2 Minimum-jerk Trajectory

close all;

% The time step of the simulation and its number of iteration
dt = 1e-3;
Nt = 6000;

% The total time and its time array
T     = dt * Nt;
t_arr = dt * (0:(Nt-1));

% First Min-jerk Trajectory with DMP
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
