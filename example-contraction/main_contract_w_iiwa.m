%% [Example Script] Contraction Theory and DMP
% [Author] Moses Chong-ook Nah
% [Date]   2023.10.08

%% [Initialization] 
fig_config( 'fontSize', 20, 'markerSize', 10 )

%%
robot = iiwa14( 'high' );
robot.init( );

anim = Animation( 'Dimension', 3, 'xLim', [-0.1, 0.8], 'yLim', [-0.4, 0.4], 'zLim', [0.0, 0.8], 'isSaveVideo', false );
anim.init( );
% anim.attachRobot( robot )  




load( 'learned_parameters/min_jerk.mat' );


data1 = data;
data4 = data;

% Load the dataset 
load( 'learned_parameters/cosine.mat'   );
data2 = data;

load( 'learned_parameters/radial.mat'   );
data3 = data;

data_whole = { data1, data2, data3, data4 };

% Set all data to origin
data1.y_arr = data1.y_arr - data1.y_arr( :, 1 );
data2.y_arr = data2.y_arr - data2.y_arr( :, 1 );
data3.y_arr = data3.y_arr - data3.y_arr( :, 1 );
data4.y_arr = data4.y_arr - data4.y_arr( :, 1 );

data1.y_arr = data1.y_arr + [ 0.5, -0.3, 0.1 ]';
data2.y_arr = data2.y_arr + data1.y_arr( :, end ) + [ 0.0,  0.1, -0.1 ]';
data3.y_arr = data3.y_arr + data2.y_arr( :, end ) + [ 0.0, -0.4,  0.1 ]';
data4.y_arr = data4.y_arr + data3.y_arr( :, end ) + [ 0.0,  0.3, -0.1 ]';

c_arr = [ 0.0000, 0.4470, 0.7410;
          0.8500, 0.3250, 0.0980;
          0.4940, 0.1840, 0.5560;
          0.4660, 0.6740, 0.1880 ];

plot3( anim.hAxes, data1.y_arr( 1, : ), data1.y_arr( 2, : ), data1.y_arr( 3, : ), 'color', c_arr( 1, : ), 'linewidth', 5 )
plot3( anim.hAxes, data2.y_arr( 1, : ), data2.y_arr( 2, : ), data2.y_arr( 3, : ), 'color', c_arr( 2, : ), 'linewidth', 5 )
plot3( anim.hAxes, data3.y_arr( 1, : ), data3.y_arr( 2, : ), data3.y_arr( 3, : ), 'color', c_arr( 3, : ), 'linewidth', 5 )
plot3( anim.hAxes, data4.y_arr( 1, : ), data4.y_arr( 2, : ), data4.y_arr( 3, : ), 'color', c_arr( 4, : ), 'linewidth', 5 )

view( 90, 0 )

data_whole{ 1 }.p0i = data1.y_arr( :, 1 );
data_whole{ 2 }.p0i = data2.y_arr( :, 1 );
data_whole{ 3 }.p0i = data3.y_arr( :, 1 );
data_whole{ 4 }.p0i = data4.y_arr( :, 1 );

data_whole{ 1 }.p0f = data1.y_arr( :, end );
data_whole{ 2 }.p0f = data2.y_arr( :, end );
data_whole{ 3 }.p0f = data3.y_arr( :, end );
data_whole{ 4 }.p0f = data4.y_arr( :, end );

%% [1B] Sequencing the Movements

% The initial time start 
ttmp0 = 1.0;
ttmp1 = data_whole{ 1 }.tau + 0.0;
ttmp2 = data_whole{ 2 }.tau + 0.0;
ttmp3 = data_whole{ 3 }.tau + 0.0;

ttmp_arr = [ ttmp0, ttmp1, ttmp2, ttmp3 ];

t0i_arr = cumsum( ttmp_arr );


% The time step of the simulation and its number of iteration
dt = 0.003;
Nt = round( sum( ttmp_arr )/dt ) + 1000;

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
k_gain = [ 5, 1, 1 ];

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
    x_arr( i, 5:7, 1 ) = data_whole{ i }.dp_func( 0 ) ;

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
    
            dx = A * x_curr( j, : )' + 1/data_whole{ j }.tau^2*[ zeros(4,1); f_arr ] + [ zeros(4,1); Az*data_whole{ j }.p0f ];
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

%% [1C] Drawing other stuffs

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

c_arr = [ 0.0000, 0.4470, 0.7410;
          0.8500, 0.3250, 0.0980;
          0.4940, 0.1840, 0.5560;
          0.4660, 0.6740, 0.1880 ];

% Start and end of the positions
mk = 100;
scatter3( anim.hAxes, x1(   1 ), y1(   1 ), z1(   1 ), mk , 'filled', 'o'     , 'markerfacecolor', c_arr( 1, : ), 'markeredgecolor', 'black' );
scatter3( anim.hAxes, x1( end ), y1( end ), z1( end ), mk , 'filled', 'square', 'markerfacecolor', c_arr( 1, : ), 'markeredgecolor', 'black' );

scatter3( anim.hAxes, x2(   1 ), y2(   1 ), z2(   1 ), mk , 'filled', 'o'     , 'markerfacecolor', c_arr( 2, : ), 'markeredgecolor', 'black' );
scatter3( anim.hAxes, x2( end ), y2( end ), z2( end ), mk , 'filled', 'square', 'markerfacecolor', c_arr( 2, : ), 'markeredgecolor', 'black' );

scatter3( anim.hAxes, x3(   1 ), y3(   1 ), z3(   1 ), mk , 'filled', 'o'     , 'markerfacecolor', c_arr( 3, : ), 'markeredgecolor', 'black' );
scatter3( anim.hAxes, x3( end ), y3( end ), z3( end ), mk , 'filled', 'square', 'markerfacecolor', c_arr( 3, : ), 'markeredgecolor', 'black' );

scatter3( anim.hAxes, x4(   1 ), y4(   1 ), z4(   1 ), mk , 'filled', 'o'     , 'markerfacecolor', c_arr( 4, : ), 'markeredgecolor', 'black' );
scatter3( anim.hAxes, x4( end ), y4( end ), z4( end ), mk , 'filled', 'square', 'markerfacecolor', c_arr( 4, : ), 'markeredgecolor', 'black' );

% Running the Animation Loop
s1 = scatter3( anim.hAxes, x1( 1 ), y1( 1 ), z1( 1 ), 1.3*mk , 'filled', 'o' , 'markerfacecolor', c_arr( 1, : ), 'markeredgecolor', 'black' );
s2 = scatter3( anim.hAxes, x2( 1 ), y2( 1 ), z2( 1 ), 1.3*mk, 'filled', 'o' , 'markerfacecolor', c_arr( 2, : ), 'markeredgecolor', 'black' );
s3 = scatter3( anim.hAxes, x3( 1 ), y3( 1 ), z3( 1 ), 1.3*mk, 'filled', 'o' , 'markerfacecolor', c_arr( 3, : ), 'markeredgecolor', 'black' );
s4 = scatter3( anim.hAxes, x4( 1 ), y4( 1 ), z4( 1 ), 1.3*mk, 'filled', 'o' , 'markerfacecolor', c_arr( 4, : ), 'markeredgecolor', 'black' );

scc = scatter3( anim.hAxes, xc( 1 ), yc( 1 ), zc( 1 ), 3.0*mk, 'filled', 'o' , 'markerfacecolor', [ 0.5, 0.5, 0.5], 'markeredgecolor', 'black', 'markerfacealpha', 0.8 );

s_arr = { s1, s2, s3, s4 };

% v = VideoWriter( 'video.mp4','MPEG-4' );
% v.FrameRate = 30;
t1 = title( anim.hAxes, sprintf( 'Time %.3f s', 0) );

% open( v );
tmp_step = 33;
for i = 1 : tmp_step : Nt
    
    for j = 1 : length( s_arr )
        set( s_arr{ j }, 'XData', x_arr_whole( i, j ), 'YData', y_arr_whole( i, j ), 'ZData', z_arr_whole( i, j ) );
    end
    
    set( scc, 'XData', xc( i + 1 ), 'YData', yc( i + 1 ), 'ZData', zc( i + 1 ))

    drawnow 
    
    set( t1, 'string', sprintf( 'Time %.3f s', t_arr( i ) ) );
    
    tmp_frame = getframe( anim.hFig );
    % writeVideo( v,tmp_frame );
    i
end
% close( v );

%% [1C-1] Save the Data

y_arr = [ xc; yc; zc ];

dir_name  = './data/example7/';

writematrix(  y_arr - y_arr( :, 1 ), [ dir_name, 'pos_data.csv' ] ) 


%% [1D] Actual Control of iiwa with xc, yc, zc

% Update kinematics
clear pi
q_init = [   -0.4271    1.3218   -0.1504   -1.2475   -0.0905    0.8930   -0.6004]';
% q_init = [0, 28.56, 0, -87.36, -7.82, 75.56, -9.01]' * pi/180;

robot.updateKinematics( q_init );
anim.update( 0 );

% Title: simulation time
titleString = sprintf( 'Time: %2.1f sec', 0 );
mytitle = title( titleString );
set( mytitle, 'FontSize' , 15);

% Point-to-point Joint Space Impedance Controller
% Using Minimum-jerk trajectory as the virtual trajectory 

% Initial joint posture and velocity
q  = q_init;
dq = zeros( robot.nq, 1 );

% The initial end-effector position.
Hi = robot.getForwardKinematics( q );
pi = Hi( 1:3, 4 );

% Get the current SO(3) orientation matrix. 
Rinit = Hi( 1:3, 1:3 ); 

% Time step for the simulation
ns = 0;

% End-effector task-space impedances
Kp = 400 * eye( 3 );
Bp = 0.1 * Kp;

Bq = 0.5 * eye( robot.nq );

t0 = 0.1;
D  = 2.0;


kr = 3;
br = 0.1*kr;
p_ref = zeros( 3, 1 );
i = 1;
t = 0;

v = VideoWriter( 'video2.mp4','MPEG-4' );
v.FrameRate = 30;
t1 = title( anim.hAxes, sprintf( 'Time %.3f s', 0) );

anim.gPatches{8}.EdgeAlpha = 0.3;
anim.gPatches{7}.EdgeAlpha = 0.3;
anim.gPatches{8}.FaceAlpha = 0.3;
anim.gPatches{7}.FaceAlpha = 0.3;

open( v );
tmp_step = 33;

ns = 1;
while t <= T
    
    % Get the mass matrix of the Acrobot
    M = robot.getMassMatrix( q );

    % Get the Coriolis term of the robot
    C = robot.getCoriolisMatrix( q, dq );
    
    % Get the Gravity term of the robot
    G = robot.getGravityVector( q );

    % Get the Hybrid Jacobian 
    JH = robot.getHybridJacobian( q );
    

    % Get the end-effector position and velocity 
    JHp = JH( 1:3, : );     % Hybrid Jacobian, end-effector position
    JHr = JH( 4:6, : );     % Hybrid Jacobian, end-effector orientation
    dp = JHp * dq;
    
    % The initial end-effector position and orientation 
    H = robot.getForwardKinematics( q );
    p = H( 1:3, 4 );
    Rcurr = H( 1:3, 1:3 );

    % Get the angle between Rcurr and Rinit
    Rdel = ( Rcurr )^-1 * Rinit;
    theta = acos( ( trace( Rdel ) - 1 )/2 );

    % Get the axis
    if theta ~= 0
        w_mat = ( Rdel - Rdel' ) / ( 2 * sin( theta ) ); 
    else
        w_mat = zeros( 3, 3 );
    end

    w = [ -w_mat( 2,3 ), w_mat( 1,3 ), -w_mat( 1,2 ) ]';

    p_ref( 1 ) = xc( i );
    p_ref( 2 ) = yc( i );
    p_ref( 3 ) = zc( i );

    tau1 = JHp' * ( Kp * ( p_ref - p ) + Bp * ( - dp ) );
    tau2 = JHr' * ( kr * Rcurr * w * theta - br * JHr * dq );

    tau = tau1 + tau2 - Bq * dq;
    rhs = M\( -C * dq + tau ); 

%     rhs = zeros( robot.nq, 1 );
    
    % We compare the C matrices with the matrices provided from Eq.9
    % The default values of our Acrobot geometrical/inertial parameters are as follows:
    % m1 = 1, m2 = 1, lc1 = 0.5, lc2 = 0.5, l1 = 1.0, l2 = 1.0    
    % The following 
    % [REF] https://underactuated.csail.mit.edu/acrobot.html#section1
    
    [ q1, dq1 ] = func_symplecticEuler( q, dq, rhs, dt );
    q  =  q1;
    dq = dq1;
    
    robot.updateKinematics( q );
    anim.update( t );    
    
  
    if round( t / anim.FrameUpdateTime ) >= ns
        % Update the linkage plot
        robot.updateKinematics( q );
        anim.update( t );    
        
        
        for j = 1 : length( s_arr )
            set( s_arr{ j }, 'XData', x_arr_whole( i, j ), 'YData', y_arr_whole( i, j ), 'ZData', z_arr_whole( i, j ) );
        end
        
        set( scc, 'XData', xc( i ), 'YData', yc( i ), 'ZData', zc( i  ) )    

        % Set animation title
        set( mytitle, 'String', sprintf( 'Time: %2.1f sec', t ) );
        ns = ns + 1;

        tmp_frame = getframe( anim.hFig );
        writeVideo( v,tmp_frame );

    end    


    % Get the forward kinematics of the EE
    t = t + dt;                                                                
    i = i + 1;
end

anim.close( )
close( v );
