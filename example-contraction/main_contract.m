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

data_whole{ 3 }.p0i = data_whole{ 2 }.p0f + [ 1, 1, 1 ];
data_whole{ 3 }.p0f = data_whole{ 3 }.p0i + data_whole{ 3 }.p0f;

data_whole{ 4 }.p0i = data_whole{ 3 }.p0f + [ 1, 3, -1];
data_whole{ 4 }.p0f = data_whole{ 4 }.p0i + data_whole{ 4 }.p0f;

%% [1B] Sequencing the Movements

close all;

% The initial time start 
ttmp0 = 1.0;
ttmp1 = data_whole{ 1 }.tau + 0.0;
ttmp2 = data_whole{ 2 }.tau + 0.0;
ttmp3 = data_whole{ 3 }.tau + 0.0;

ttmp_arr = [ ttmp0, ttmp1, ttmp2, ttmp3 ];

t0i_arr = cumsum( ttmp_arr );


% The time step of the simulation and its number of iteration
dt = 1e-3;
Nt = round( sum( ttmp_arr )/dt ) + 5000;

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

x_coupled_curr = zeros( 7, 1 );
x_coupled  = zeros( 7, Nt );
dx_coupled = zeros( 7, Nt );

% The State-array 
Az_arr = zeros( 4, 3, 3 );
Bz_arr = zeros( 4, 3, 3 );
A_arr  = zeros( 4, 7, 7 );

% k_gains
k_gain = [ 5, 5, 1 ];

% Timing for activation
% ttmp_arr = [ ttmp0, ttmp1, ttmp2, ttmp3 ];
t_start_arr = [ t0i_arr( 1 ) + 0.8*data_whole{ 1 }.tau, ...
                t0i_arr( 2 ) + 0.8*data_whole{ 2 }.tau, ...
                t0i_arr( 3 ) + 0.8*data_whole{ 3 }.tau ];

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

x_curr = x0_arr;
t = 0;
act_arr = zeros( 1, 3 );

% Forward Integration
for i = 0 : (Nt-1)

    % Integrating over the movement Separately
    for j = 1 : 4
        
        if t <= t0i_arr( j )
            x_arr( j, :, i + 1 )  = x0_arr( j, : );
            x_curr( j, : )        = x0_arr( j, : );
            
            if j == 1
                x_coupled( :, i + 1 ) = x0_arr( j, : );   
                x_coupled_curr        = x0_arr( j, : )';   
            end
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
    
            dx = A * x_curr( j, : )' + 1/data_whole{ j }.tau^2*[ zeros(4,1); f_arr ] + [ zeros(4,1); Az*data_whole{ j }.p0f' ];
            x_curr( j, : ) = x_curr( j, : ) + dx' * dt;
            x_arr( j, :, i+1 ) = x_curr( j, : );    

            % For the first movement, add the dynamics of the coupling
                % For the transition time 
                if     t_start_arr( 1 ) <= t && t_start_arr( 2 ) >= t
                      act_arr( 1 ) = act_arr( 1 ) + 0.001;
                      if act_arr( 1 ) >= 1 
                        act_arr( 1 ) = 1;
                      end
                      dx_coupled = dx - act_arr( 1 ) * k_gain( 1 ) * ( x_coupled_curr - x_curr( 2, : )' );

                elseif t_start_arr( 2 ) <= t && t_start_arr( 3 ) >= t
                      act_arr( 2 ) = act_arr( 2 ) + 0.001;
                      if act_arr( 2 ) >= 1 
                        act_arr( 2 ) = 1;
                      end                    
                      dx_coupled = dx - act_arr( 2 ) * k_gain( 2 ) * ( x_coupled_curr - x_curr( 3, : )' );

                elseif t <= t_start_arr( 1 )
                      dx_coupled = dx;
                else
                      act_arr( 3 ) = act_arr( 3 ) + 0.001;
                      if act_arr( 3 ) >= 1 
                        act_arr( 3 ) = 1;
                      end 
                      dx_coupled = dx - act_arr( 3 ) * k_gain( 3 ) * ( x_coupled_curr - x_curr( 4, : )' );
                end
                
                x_coupled_curr = x_coupled_curr + dx_coupled * dt;
                x_coupled( :, i+1 ) = x_coupled_curr;


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

xc = squeeze( x_coupled( 2, : ) );
yc = squeeze( x_coupled( 3, : ) );
zc = squeeze( x_coupled( 4, : ) );

x_arr_whole = [ x1, x2, x3, x4 ];
y_arr_whole = [ y1, y2, y3, y4 ];
z_arr_whole = [ z1, z2, z3, z4 ];

f = figure( ); a = axes( 'parent', f );
hold on;

c_arr = [ 0.0000, 0.4470, 0.7410;
          0.8500, 0.3250, 0.0980;
          0.4940, 0.1840, 0.5560;
          0.4660, 0.6740, 0.1880 ];

plot3( x1, y1, z1, 'color', c_arr( 1, : ), 'linewidth', 3 )
plot3( x2, y2, z2, 'color', c_arr( 2, : ), 'linewidth', 3 )
plot3( x3, y3, z3, 'color', c_arr( 3, : ), 'linewidth', 3 )
plot3( x4, y4, z4, 'color', c_arr( 4, : ), 'linewidth', 3 )



% Start and end of the positions
scatter3( x1(   1 ), y1(   1 ), z1(   1 ), 200, 'filled', 'o'     , 'markerfacecolor', c_arr( 1, : ), 'markeredgecolor', 'black' );
scatter3( x1( end ), y1( end ), z1( end ), 200, 'filled', 'square', 'markerfacecolor', c_arr( 1, : ), 'markeredgecolor', 'black' );

scatter3( x2(   1 ), y2(   1 ), z2(   1 ), 200, 'filled', 'o'     , 'markerfacecolor', c_arr( 2, : ), 'markeredgecolor', 'black' );
scatter3( x2( end ), y2( end ), z2( end ), 200, 'filled', 'square', 'markerfacecolor', c_arr( 2, : ), 'markeredgecolor', 'black' );

scatter3( x3(   1 ), y3(   1 ), z3(   1 ), 200, 'filled', 'o'     , 'markerfacecolor', c_arr( 3, : ), 'markeredgecolor', 'black' );
scatter3( x3( end ), y3( end ), z3( end ), 200, 'filled', 'square', 'markerfacecolor', c_arr( 3, : ), 'markeredgecolor', 'black' );

scatter3( x4(   1 ), y4(   1 ), z4(   1 ), 200, 'filled', 'o'     , 'markerfacecolor', c_arr( 4, : ), 'markeredgecolor', 'black' );
scatter3( x4( end ), y4( end ), z4( end ), 200, 'filled', 'square', 'markerfacecolor', c_arr( 4, : ), 'markeredgecolor', 'black' );

lw = 10;
tmp_s = [ 0, 0, 0];
set( a, 'view', [65.1000, 11.7203], 'xlim', [tmp_s( 1 ), tmp_s( 1 )+lw], ...
                                    'ylim', [tmp_s( 2 ), tmp_s( 2 )+lw], ...
                                    'zlim', [tmp_s( 3 ), tmp_s( 3 )+lw] )
set( a, 'fontsize', 30 )

% Running the Animation Loop
s1 = scatter3( x1( 1 ), y1( 1 ), z1( 1 ), 400, 'filled', 'o' , 'markerfacecolor', c_arr( 1, : ), 'markeredgecolor', 'black' );
s2 = scatter3( x2( 1 ), y2( 1 ), z2( 1 ), 400, 'filled', 'o' , 'markerfacecolor', c_arr( 2, : ), 'markeredgecolor', 'black' );
s3 = scatter3( x3( 1 ), y3( 1 ), z3( 1 ), 400, 'filled', 'o' , 'markerfacecolor', c_arr( 3, : ), 'markeredgecolor', 'black' );
s4 = scatter3( x4( 1 ), y4( 1 ), z4( 1 ), 400, 'filled', 'o' , 'markerfacecolor', c_arr( 4, : ), 'markeredgecolor', 'black' );

scc = scatter3( xc( 1 ), yc( 1 ), zc( 1 ), 500, 'filled', 'o' , 'markerfacecolor', [ 0.5, 0.5, 0.5], 'markeredgecolor', 'black', 'markerfacealpha', 0.5 );

s_arr = { s1, s2, s3, s4 };

v = VideoWriter( 'video.mp4','MPEG-4' );
v.FrameRate = 30;
t1 = title( sprintf( 'Time %.3f s', 0) );

open( v );
tmp_step = 33;
for i = 1 : tmp_step : Nt
    
    for j = 1 : length( s_arr )
        set( s_arr{ j }, 'XData', x_arr_whole( i, j ), 'YData', y_arr_whole( i, j ), 'ZData', z_arr_whole( i, j ) );
    end
    
    set( scc, 'XData', xc( i + 1 ), 'YData', yc( i + 1 ), 'ZData', zc( i + 1 ))

    drawnow 
    
    set( t1, 'string', sprintf( 'Time %.3f s', t_arr( i ) ) );
    
    tmp_frame = getframe( f );
    writeVideo( v,tmp_frame );
    i
end
close( v );