% [Title]     Generating Images for the revised manuscript of IJRR
% [Author]    Moses Chong-ook Nah
% [Email]     mosesnah@mit.edu
% [Update]    At 2024.01.11

%% (--) Initialization
clear; close all; clc;
fig_config( 'fontSize', 20, 'markerSize', 10 )

%% =======================================================
%% (1-) Data Processing, for Drawing M
%% (--) (1A) Checking Data

% Under example4 directory under data
raw_data = parse_txt( 'data/example4/iiwa_example_pos.txt', 0 );

% Read the time (with offset) and joint-trajectory data
t_arr_demo = raw_data( :, 1 )' - raw_data( 1, 1 );
q_arr_demo = raw_data( :, 2:8 )';

% Call the robot for Forward Kinematics
robot = iiwa14( 'high' );
robot.init( );

% The number of sample points P, which we denote as Np
[ ~, Np ] = size( q_arr_demo );

% End-effect Forward Kinematics, both position and orientation
% We will learn the position, but also saving the orientation for backup
p_arr_raw = zeros( 3, Np );
R_arr_raw = zeros( 3, 3, Np );

% Calculating the position/orientation with Forward Kinematics
for i = 1 : Np
    H = robot.getForwardKinematics( q_arr_demo( :, i ) );
    p_arr_raw( :, i )    = H( 1:3,   4 );
    R_arr_raw( :, :, i ) = H( 1:3, 1:3 ); 
end

% Minus the offset for initial dataset zero 
p_arr_raw = p_arr_raw - p_arr_raw( :, 1 );

% Simple time diff for comparison
dp_arr_raw  = data_diff(  p_arr_raw );
ddp_arr_raw = data_diff( dp_arr_raw );

% Filtering out the Data
  p_arr_filt = zeros( size( p_arr_raw ), 'like', p_arr_raw );
 dp_arr_filt = zeros( size( p_arr_raw ), 'like', p_arr_raw );
ddp_arr_filt = zeros( size( p_arr_raw ), 'like', p_arr_raw );

% Filtered position, velocity and acceleration
for i = 1 : 3
    % Gaussian Filtering of Position 
     p_arr_filt( i, : ) = smoothdata( p_arr_raw( i, : ), "gaussian", 50 );

    % Diff and Gaussian Filter of Velocity  
   tmp_vel = data_diff( p_arr_filt );
    dp_arr_filt( i, : ) = smoothdata( tmp_vel( i, : ), "gaussian", 50 );

    % Diff and Gaussian Filter of Acceleration
   tmp_acc = data_diff( dp_arr_filt );
   ddp_arr_filt( i, : ) = smoothdata( tmp_acc( i, : ), "gaussian", 50 );   
end


% Plotting the Data
tmp_plot = cell( 1, 6 );
tmp_plot{ 1 } =   p_arr_raw; tmp_plot{ 2 } =   p_arr_filt;
tmp_plot{ 3 } =  dp_arr_raw; tmp_plot{ 4 } =  dp_arr_filt;
tmp_plot{ 5 } = ddp_arr_raw; tmp_plot{ 6 } = ddp_arr_filt;

tmp_label = [ "Pos.", "Vel.", "Acc." ];

% XYZ Color
c_arr = [0.0000, 0.4470, 0.7410;
         0.8500, 0.3250, 0.0980;
         0.9290, 0.6940, 0.1250];

% Row
for i = 1:3

    % Column
    for j = 1 :2
        subplot( 3, 2, 2*(i-1)+j )
        hold on
        
        % Iterating over XYZ
        for k = 1 : 3
            plot( t_arr_demo, tmp_plot{ 2*(i-1)+j }( k, : ), 'linewidth', 3, 'color', c_arr( k, : ) )
        end
        
        % Ignore the first and second row's time 
        if i == 1 || i == 2
            set( gca, 'xticklabel', {} )
        end
        
        % Add time label
        if i == 3
           xlabel( '$t (s)$', 'fontsize', 40 )
        end
        set( gca, 'xlim', [ 0, max( t_arr_demo )+0.001 ], 'fontsize', 30 )
        
        if j == 1
           ylabel( tmp_label( i ), 'fontsize', 40 )
           tmp = get( gca, 'ylim' );
        else
           set( gca, 'ylim', tmp ) 
        end
    end
    
end

mySaveFig( gcf, 'revision_MIM_pos' )


%% (--) (1B) Imaging for learning Drawing M
% For this, iiwa14 is required.
close all; clear all;
% Under example4 directory under data
raw_data = parse_txt( 'data/example4/iiwa_example_pos.txt', 0 );

% Read the time (with offset) and joint-trajectory data
t_arr_demo = raw_data( :, 1 )' - raw_data( 1, 1 );
q_arr_demo = raw_data( :, 2:8 )';

% Call the robot for Forward Kinematics
robot = iiwa14( 'high' );
robot.init( );

% The number of sample points P, which we denote as Np
[ ~, Np ] = size( q_arr_demo );

% Need animation for Visualization
anim = Animation( 'Dimension', 3, 'xLim', [-0.9,0.9], 'yLim', [-0.9,0.9], 'zLim', [0,1.2] );
anim.init( );

view( [ 90, 0 ] )

% End-effect Forward Kinematics, both position and orientation
% We will learn the position, but also saving the orientation for backup
p_arr_raw = zeros( 3, Np );
R_arr_raw = zeros( 3, 3, Np );

% Calculating the position/orientation with Forward Kinematics
for i = 1 : Np
    H = robot.getForwardKinematics( q_arr_demo( :, i ) );
    p_arr_raw( :, i )    = H( 1:3,   4 );
    R_arr_raw( :, :, i ) = H( 1:3, 1:3 ); 
end

set( anim.hAxes, 'visible', 'off' )
% Overlap the M drawing

y_off    = [-0.5, 0.0, 0.5];
x = 0.1;

for i = 1 :length( y_off )
    y = y_off( i );

    plot3( anim.hAxes, x+ p_arr_raw( 1, : ), y + p_arr_raw( 2, : ), p_arr_raw( 3, : ), 'linewidth', 10, 'color',[0.4940 0.1840 0.5560]	, 'linestyle', ':' )
    scatter3( anim.hAxes,x+ p_arr_raw( 1,   1 ), y + p_arr_raw( 2,   1 ), p_arr_raw( 3,   1 ), 300, 'filled', 'd', 'markerfacecolor', 'white', 'markeredgecolor', 'k', 'linewidth', 5)
    scatter3( anim.hAxes,x+ p_arr_raw( 1, end ), y + p_arr_raw( 2, end ), p_arr_raw( 3, end ), 300, 'filled', 'o', 'markerfacecolor', 'white', 'markeredgecolor', 'k', 'linewidth', 5)
    
end


% Setting some specific Time lines that we want to plot 
N_moment = [ 1, 800, Np ];
b_alpha  = [ 0.5, 0.5, 0.5];

for i = 1 : length( N_moment )
    
    % Call the robot for Forward Kinematics
    robot = iiwa14( 'high' );
    robot.init( );
    
    tmp = eye( 4 );
    tmp( 2, 4 ) = y_off( i );
    anim.Robots{ i }.H_base = tmp;    
    robot.H_base = tmp;
    
    anim.attachRobot( robot )      
    
    y = y_off( i );

    N = N_moment( i );
    b = b_alpha( i );
    q = q_arr_demo( :, N );
    
    anim.Robots{ i }.updateKinematics( q );
    for j = 1 : 8
        if ismember( j , 1:6 )
            set( anim.gPatches{ i }{ j }, 'FaceAlpha', 0.2 )
        else
            set( anim.gPatches{ i }{ j }, 'FaceAlpha', b )
        end
    end
    
    scatter3( anim.hAxes, x+ p_arr_raw( 1, N ), y + p_arr_raw( 2, N ), p_arr_raw( 3, N ), 300, 'filled', 'o', 'markerfacecolor', 'white', 'markeredgecolor', 'k', 'linewidth', 5)
    
end

anim.update( 0 );
mySaveFig( gcf, 'revision_MIM_pos_KUKA_iiwa14' )


%% (--) (1C) Actual Result for this Input
close all; clc;

% The actual result
% Under example4 directory under data
raw_data  = parse_txt( 'data/example4/draw_M_actual_data1.txt', 0 );

% Read the time (with offset) and joint-trajectory data
t_arr_demo  = raw_data( :, 1 )' - raw_data( 1, 1 );
q_arr_demo  = raw_data( :, 2:8 )';
p_arr_demo  = raw_data( :, 9:11 )';
p0_arr_demo = raw_data( :, 13:end )';

% Call the robot for Forward Kinematics
robot = iiwa14( 'high' );
robot.init( );

% The number of sample points P, which we denote as Np
[ ~, Np ] = size( q_arr_demo );

% Need animation for Visualization
anim = Animation( 'Dimension', 3, 'xLim', [-0.9,0.9], 'yLim', [-0.9,0.9], 'zLim', [0,1.2] );
anim.init( );
set( anim.hAxes, 'visible', 'off' )


view( [ 90, 0 ] )


% Setting some specific Time lines that we want to plot 
N_moment = [ 1, 1800,  Np ];
b_alpha  = [ 0.6, 0.6, 0.6];

y_off    = [-0.5, 0.0, 0.5];
x = 0.1;

for i = 1:length( y_off )
    y = y_off( i );

    plot3( anim.hAxes, x+ p_arr_demo( 1, : ), y + p_arr_demo( 2, : ), p_arr_demo( 3, : ), 'linewidth', 10, 'color', 'k', 'linestyle', '-' )
    plot3( anim.hAxes, x+ p0_arr_demo( 1, : ), y + p0_arr_demo( 2, : ), p0_arr_demo( 3, : ), 'linewidth', 7, 'color',[0.4940 0.1840 0.5560]	, 'linestyle', ':' )

    scatter3( anim.hAxes,x+ p_arr_demo( 1,   1 ), y + p_arr_demo( 2,   1 ), p_arr_demo( 3,   1 ), 300, 'filled', 'd', 'markerfacecolor', 'white', 'markeredgecolor', 'k', 'linewidth', 5)
    scatter3( anim.hAxes,x+ p_arr_demo( 1, end ), y + p_arr_demo( 2, end ), p_arr_demo( 3, end ), 300, 'filled', 'o', 'markerfacecolor', 'white', 'markeredgecolor', 'k', 'linewidth', 5)
    
end

for i = 1 : length( N_moment )
    
    % Call the robot for Forward Kinematics
    robot = iiwa14( 'high' );
    robot.init( );
    
    tmp = eye( 4 );
    tmp( 2, 4 ) = y_off( i );
    anim.Robots{ i }.H_base = tmp;    
    robot.H_base = tmp;
    y = y_off( i );
       
    
    anim.attachRobot( robot )  
    
    N = N_moment( i );
    b = b_alpha( i );
    q = q_arr_demo( :, N );
   
    anim.Robots{ i }.updateKinematics( q );
    for j = 1 : 8
        if ismember( j , 1:6 )
            set( anim.gPatches{ i }{ j }, 'FaceAlpha', 0.2 )
        else
            set( anim.gPatches{ i }{ j }, 'FaceAlpha', b )
        end
    end

    scatter3( anim.hAxes, x+ p_arr_demo( 1, N ), y + p_arr_demo( 2, N ), p_arr_demo( 3, N ), 300, 'filled', 'o', 'markerfacecolor', 'white', 'markeredgecolor', 'k', 'linewidth', 5)
        
end
anim.update( 0 )
mySaveFig( gcf, 'revision_MIM_pos_result' )


%% =======================================================
%% (2-) Data Processing, for Lifting Up and Down
%% (--) (2A) Plotting the Quaternions 
close all;clc;
% Under example4 directory under data
raw_data = parse_txt( 'data/example3/iiwa_example_orient.txt', 0 );


% Read the time (with offset) and joint-trajectory data
t_arr_demo = raw_data( :, 1 )' - raw_data( 1, 1 );
q_arr_demo = raw_data( :, 2:8 )';

% Call the robot for Forward Kinematics
robot = iiwa14( 'high' );
robot.init( );

% The number of sample points P, which we denote as Np
[ ~, Np ] = size( q_arr_demo );

% Get the raw quaternion matrix 
R_arr_raw    = zeros( 3, 3, Np );
quat_arr_raw = zeros( 4, Np );

for i = 1 : Np
    H = robot.getForwardKinematics( q_arr_demo( :, i ) );
    R_arr_raw( :, :, i ) = H( 1:3, 1:3 ); 
    quat_arr_raw( :, i ) = rotm2quat( H( 1:3, 1:3 ) );
end

f = figure( ); a = axes( 'parent', f );
plot( a, t_arr_demo, quat_arr_raw, 'linewidth', 6, 'color', 'k' )
set( a, 'fontsize', 40, 'xlim', [ 0, max( t_arr_demo )], 'xtick', 0:1:max( t_arr_demo ), 'ytick', -0.8:0.4:0.4001 )
xlabel( '$t (s)$', 'fontsize', 40 )
ylabel( '$\vec{\mathbf{q}}(t)$', 'fontsize', 40 )

mySaveFig( gcf, 'revision_MIM_get_quat' )

%% (--) (2B) Checking Data, Error Data
close all;clc;

% Under example4 directory under data
raw_data = parse_txt( 'data/example3/iiwa_example_orient.txt', 0 );

% Read the time (with offset) and joint-trajectory data
t_arr_demo = raw_data( :, 1 )' - raw_data( 1, 1 );
q_arr_demo = raw_data( :, 2:8 )';

% Call the robot for Forward Kinematics
robot = iiwa14( 'high' );
robot.init( );

% The number of sample points P, which we denote as Np
[ ~, Np ] = size( q_arr_demo );

% Get the raw quaternion matrix 
R_arr_raw    = zeros( 3, 3, Np );
quat_arr_raw = zeros( 4, Np );

for i = 1 : Np
    H = robot.getForwardKinematics( q_arr_demo( :, i ) );
    R_arr_raw( :, :, i ) = H( 1:3, 1:3 ); 
    quat_arr_raw( :, i ) = rotm2quat( H( 1:3, 1:3 ) );
end

quat_goal = quat_arr_raw( :, end );
error_arr_raw  = zeros( 3, Np );
error_arr_filt = zeros( 3, Np );
for i = 1 : Np
    error_arr_raw( :, i ) = get_quat_error( quat_arr_raw( :, i ), quat_goal  );
end

for i = 1 : 3
    error_arr_filt( i, : ) = smoothdata( error_arr_raw( i, : ), "gaussian", 50 );
end

% Diff and Gaussian Filter of Velocity  
dq_arr_demo = zeros( 7, Np );
for i = 1 : 7
    dq_arr_demo( i, : ) = data_diff( q_arr_demo( i, : ) );
end

% Plot the another derivative for the error array 
ws_arr    = zeros( 3, Np );
dquat_arr = zeros( 4, Np );
for i = 1 : Np
    JS  = robot.getSpatialJacobian( q_arr_demo( :, i ) );
    JSr = JS( 4:6, : );
    ws_arr( :, i ) = JSr * dq_arr_demo( :, i );
    dquat_arr( :, i ) = 1/2 * quat_mul( R3_to_quat( ws_arr( :, i ) ), quat_arr_raw( :, i  )' );
end

derror_arr_raw  = zeros( 3, Np );
derror_arr_filt = zeros( 3, Np );

for i = 1 : Np
    derror_arr_raw( :, i ) = dLogQuat( quat_arr_raw( :, i ), dquat_arr( :, i ) );
end

for i = 1 : 3
    derror_arr_filt( i, : ) = smoothdata( derror_arr_raw( i, : ), "gaussian", 50 );
end

dderror_arr_raw  = zeros( 3, Np );
dderror_arr_filt = zeros( 3, Np );

for i = 1 : 3
    dderror_arr_raw(  i, : ) = data_diff(  derror_arr_raw( i, : ) );
    dderror_arr_filt( i, : ) = smoothdata( dderror_arr_raw( i, : ), "gaussian", 50 );
end

% Plotting the Data
tmp_plot = cell( 1, 6 );
tmp_plot{ 1 } =   error_arr_raw; tmp_plot{ 2 } =   error_arr_filt;
tmp_plot{ 3 } =  derror_arr_raw; tmp_plot{ 4 } =  derror_arr_filt;
tmp_plot{ 5 } = dderror_arr_raw; tmp_plot{ 6 } = dderror_arr_filt;

c_arr = [0.0000, 0.4470, 0.7410;
         0.8500, 0.3250, 0.0980;
         0.9290, 0.6940, 0.1250];
     
ylim_arr = [ -0.5, 1.5; -4e-3, 4e-3; -6e-5, 6e-5];
ytick_arr =  cell( 1, 3 );
ytick_arr{ 1 } = [-0.5, 0.5, 1.5];
ytick_arr{ 2 } = [-4e-3, 0, 4e-3];
ytick_arr{ 3 } = [-5e-5, 0, 5e-5];

tmp_label = [ "$\mathbf{e}(t)$", "$\dot{\mathbf{e}}(t)$", "${\ddot{\mathbf{e}}(t)}$" ];
% Row
for i = 1:3

    % Column
    for j = 1 :2
        subplot( 3, 2, 2*(i-1)+j )
        hold on
        
        % Iterating over XYZ
        for k = 1 : 3
            plot( t_arr_demo, tmp_plot{ 2*(i-1)+j }( k, : ), 'linewidth', 3, 'color', c_arr( k, : ) )
        end
        
        % Ignore the first and second row's time 
        if i == 1 || i == 2
            set( gca, 'xticklabel', {} )
        end
        
        % Add time label
        if i == 3
           xlabel( '$t (s)$', 'fontsize', 40 )
        end
        set( gca, 'xlim', [ 0, max( t_arr_demo )+0.001 ], 'fontsize', 30 )
        set( gca, 'ylim', ylim_arr( i, : ), 'ytick', ytick_arr{ i } )
   
        if j == 1
           ylabel( tmp_label( i ), 'fontsize', 40 )
        end        
        
    end
    
end

mySaveFig( gcf, 'revision_MIM_orient' )

%% (--) (2C) KUKA iiwa14 Learning Process, Lifting up and Down

% For this, iiwa14 is required.
close all; clear all;
% Under example4 directory under data
raw_data = parse_txt( 'data/example3/iiwa_example_orient.txt' , 0 );

% Read the time (with offset) and joint-trajectory data
t_arr_demo = raw_data( :, 1 )' - raw_data( 1, 1 );
q_arr_demo = raw_data( :, 2:8 )';

% Call the robot for Forward Kinematics
robot = iiwa14( 'high' );
robot.init( );

% The number of sample points P, which we denote as Np
[ ~, Np ] = size( q_arr_demo );

% Need animation for Visualization
anim = Animation( 'Dimension', 3, 'xLim', [-0.9,0.9], 'yLim', [-0.9,0.9], 'zLim', [0,1.2] );
anim.init( );


view( [ 90, 0 ] )

% End-effect Forward Kinematics, both position and orientation
% We will learn the position, but also saving the orientation for backup
p_arr_raw = zeros( 3, Np );
R_arr_raw = zeros( 3, 3, Np );

% Calculating the position/orientation with Forward Kinematics
for i = 1 : Np
    H = robot.getForwardKinematics( q_arr_demo( :, i ) );
    p_arr_raw( :, i )    = H( 1:3,   4 );
    R_arr_raw( :, :, i ) = H( 1:3, 1:3 ); 
end

set( anim.hAxes, 'visible', 'off' )
% Overlap the M drawing

% Setting some specific Time lines that we want to plot 
N_moment = [ 1, 250, 500, 750, Np ];
b_alpha  = [ 0.5, 0.5, 0.5, 0.5,0.5];
y_off    = [-0.6, -0.3, 0.0, 0.3, 0.6];
x = 0.1;

for i = 1 :length( y_off )
    y = y_off( i );

    plot3( anim.hAxes, x+ p_arr_raw( 1, : ), y + p_arr_raw( 2, : ), p_arr_raw( 3, : ), 'linewidth', 5, 'color', 'k', 'linestyle', ':' )
    scatter3( anim.hAxes,x+ p_arr_raw( 1,   1 ), y + p_arr_raw( 2,   1 ), p_arr_raw( 3,   1 ), 300, 'filled', 'd', 'markerfacecolor', 'white', 'markeredgecolor', 'k', 'linewidth', 5)
    scatter3( anim.hAxes,x+ p_arr_raw( 1, end ), y + p_arr_raw( 2, end ), p_arr_raw( 3, end ), 300, 'filled', 'o', 'markerfacecolor', 'white', 'markeredgecolor', 'k', 'linewidth', 5)
    
end

scl = 0.1;
for i = 1 : length( N_moment )
    
    % Call the robot for Forward Kinematics
    robot = iiwa14( 'high' );
    robot.init( );
    
    tmp = eye( 4 );
    tmp( 2, 4 ) = y_off( i );
    anim.Robots{ i }.H_base = tmp;    
    robot.H_base = tmp;
    
    anim.attachRobot( robot )  
    
    N = N_moment( i );
    b = b_alpha( i );
    q = q_arr_demo( :, N );
    
    anim.Robots{ i }.updateKinematics( q );

    for j = 1 : 8
        if ismember( j , 1:8 )
            set( anim.gPatches{ i }{ j }, 'FaceAlpha', 0.2 )
        else
            set( anim.gPatches{ i }{ j }, 'FaceAlpha', b )
        end
    end
    
    y = y_off( i );
    scatter3( anim.hAxes,x+ p_arr_raw( 1, N), y + p_arr_raw( 2,N ), p_arr_raw( 3,N ), 100, 'filled', 'o', 'markerfacecolor', 'white', 'markeredgecolor', [0,0,0], 'linewidth', 5)
    
    xtmp = x+ p_arr_raw( 1, N );
    ytmp = y+ p_arr_raw( 2, N );
    ztmp =    p_arr_raw( 3, N );
    
    r1 = scl * R_arr_raw( :, 1, N );
    r2 = scl * R_arr_raw( :, 2, N );
    r3 = scl * R_arr_raw( :, 3, N );
    quiver3( anim.hAxes, xtmp, ytmp, ztmp, r1( 1 ), r1( 2 ), r1( 3 ), 'linewidth', 9, 'color', 'r' );
    quiver3( anim.hAxes, xtmp, ytmp, ztmp, r2( 1 ), r2( 2 ), r2( 3 ), 'linewidth', 9, 'color', 'g' );
    quiver3( anim.hAxes, xtmp, ytmp, ztmp, r3( 1 ), r3( 2 ), r3( 3 ), 'linewidth', 9, 'color', 'b' );
    
  
end

anim.update( 0 );
mySaveFig( gcf, 'revision_MIM_orient_KUKA_iiwa14' )

%% (--) (2D) KUKA iiwa14 Actual Control


% For this, iiwa14 is required.
close all; clear all;
% Under example4 directory under data
raw_data = parse_txt( [ 'data/example8/only_orient.txt' ], 0 );

% Read the time (with offset) and joint-trajectory data
t_arr_demo = raw_data( :, 1 )' - raw_data( 1, 1 );
q_arr_demo = raw_data( :, 2:8 )';

% Call the robot for Forward Kinematics
robot = iiwa14( 'high' );
robot.init( );

% The number of sample points P, which we denote as Np
[ ~, Np ] = size( q_arr_demo );

% Need animation for Visualization
anim = Animation( 'Dimension', 3, 'xLim', [-0.9,0.9], 'yLim', [-0.9,0.9], 'zLim', [0,1.2] );
anim.init( );


view( [ 90, 0 ] )

% End-effect Forward Kinematics, both position and orientation
% We will learn the position, but also saving the orientation for backup
p_arr_raw = zeros( 3, Np );
R_arr_raw = zeros( 3, 3, Np );

% Calculating the position/orientation with Forward Kinematics
for i = 1 : Np
    H = robot.getForwardKinematics( q_arr_demo( :, i ) );
    p_arr_raw( :, i )    = H( 1:3,   4 );
    R_arr_raw( :, :, i ) = H( 1:3, 1:3 ); 
end

set( anim.hAxes, 'visible', 'off' )
% Overlap the M drawing

% Setting some specific Time lines that we want to plot 
N_moment = [ 1, 1350, 1600, 1900, Np ];
b_alpha  = [ 0.5, 0.5, 0.5, 0.5,0.5];
y_off    = [-0.6, -0.3, 0.0, 0.3, 0.6];
x = 0.1;

for i = 1 :length( y_off )
    y = y_off( i );

    plot3( anim.hAxes, x+ p_arr_raw( 1, : ), y + p_arr_raw( 2, : ), p_arr_raw( 3, : ), 'linewidth', 5, 'color', 'k', 'linestyle', ':' )
    
end

scl = 0.1;
for i = 1 : length( N_moment )
    
    % Call the robot for Forward Kinematics
    robot = iiwa14( 'high' );
    robot.init( );
    
    tmp = eye( 4 );
    tmp( 2, 4 ) = y_off( i );
    anim.Robots{ i }.H_base = tmp;    
    robot.H_base = tmp;
    
    anim.attachRobot( robot )  
    
    N = N_moment( i );
    b = b_alpha( i );
    q = q_arr_demo( :, N );
    
    anim.Robots{ i }.updateKinematics( q );

    for j = 1 : 8
        if ismember( j , 1:8 )
            set( anim.gPatches{ i }{ j }, 'FaceAlpha', 0.2 )
        else
            set( anim.gPatches{ i }{ j }, 'FaceAlpha', b )
        end
    end
    
    y = y_off( i );
    scatter3( anim.hAxes,x+ p_arr_raw( 1, N), y + p_arr_raw( 2,N ), p_arr_raw( 3,N ), 100, 'filled', 'o', 'markerfacecolor', 'white', 'markeredgecolor', [0,0,0], 'linewidth', 5)
    
    xtmp = x+ p_arr_raw( 1, N );
    ytmp = y+ p_arr_raw( 2, N );
    ztmp =    p_arr_raw( 3, N );
    
    r1 = scl * R_arr_raw( :, 1, N );
    r2 = scl * R_arr_raw( :, 2, N );
    r3 = scl * R_arr_raw( :, 3, N );
    quiver3( anim.hAxes, xtmp, ytmp, ztmp, r1( 1 ), r1( 2 ), r1( 3 ), 'linewidth', 9, 'color', 'r' );
    quiver3( anim.hAxes, xtmp, ytmp, ztmp, r2( 1 ), r2( 2 ), r2( 3 ), 'linewidth', 9, 'color', 'g' );
    quiver3( anim.hAxes, xtmp, ytmp, ztmp, r3( 1 ), r3( 2 ), r3( 3 ), 'linewidth', 9, 'color', 'b' );
    
  
end

anim.update( 0 );
mySaveFig( gcf, 'revision_MIM_orient_actual_KUKA_iiwa14' )
%% =======================================================
%% (3-) Modular Imitation Learning
%% (--) (3-) KUKA, modular Result

% For this, iiwa14 is required.
close all; clear all;
% Under example4 directory under data
raw_data = parse_txt( 'data/example8/modular_input.txt', 0 );

% Read the time (with offset) and joint-trajectory data
t_arr_demo = raw_data( :, 1 )' - raw_data( 1, 1 );
q_arr_demo = raw_data( :, 2:8 )';
p_arr_demo  = raw_data( :, 9:11 )';
p0_arr_demo = raw_data( :, 13:end )';

% Call the robot for Forward Kinematics
robot = iiwa14( 'high' );
robot.init( );

% The number of sample points P, which we denote as Np
[ ~, Np ] = size( q_arr_demo );

% Need animation for Visualization
anim = Animation( 'Dimension', 3, 'xLim', [-1.5,1.5], 'yLim', [-1.5,1.5], 'zLim', [0,1.2] );
anim.init( );

view( [ 90, 0 ] )

% End-effect Forward Kinematics, both position and orientation
% We will learn the position, but also saving the orientation for backup
p_arr_raw = zeros( 3, Np );
R_arr_raw = zeros( 3, 3, Np );

% Calculating the position/orientation with Forward Kinematics
for i = 1 : Np
    H = robot.getForwardKinematics( q_arr_demo( :, i ) );
    p_arr_raw( :, i )    = H( 1:3,   4 );
    R_arr_raw( :, :, i ) = H( 1:3, 1:3 ); 
end

set( anim.hAxes, 'visible', 'off' )
% Overlap the M drawing

% Setting some specific Time lines that we want to plot 
N_moment = [ 1, 2800, 3450, 4000, Np ];
b_alpha  = [ 0.6, 0.6, 0.6, 0.6, 0.6];
y_off    = [-1.0, -0.5, 0.0, 0.5, 1.0];
x = 0.1;

for i = 1:length( y_off )
    y = y_off( i );

    plot3( anim.hAxes, x+ p_arr_demo( 1, : ), y + p_arr_demo( 2, : ), p_arr_demo( 3, : ), 'linewidth', 6, 'color', 'k', 'linestyle', '-' )
    plot3( anim.hAxes, x+ p0_arr_demo( 1, : ), y + p0_arr_demo( 2, : ), p0_arr_demo( 3, : ), 'linewidth', 4, 'color',[0.4940 0.1840 0.5560]	, 'linestyle', ':' )

    scatter3( anim.hAxes,x+ p_arr_demo( 1,   1 ), y + p_arr_demo( 2,   1 ), p_arr_demo( 3,   1 ), 300, 'filled', 'd', 'markerfacecolor', 'white', 'markeredgecolor', 'k', 'linewidth', 5)
    scatter3( anim.hAxes,x+ p_arr_demo( 1, end ), y + p_arr_demo( 2, end ), p_arr_demo( 3, end ), 300, 'filled', 'o', 'markerfacecolor', 'white', 'markeredgecolor', 'k', 'linewidth', 5)
    
end


scl = 0.1;
for i = 1 : length( N_moment )
    
    
    % Call the robot for Forward Kinematics
    robot = iiwa14( 'high' );
    robot.init( );
    
    tmp = eye( 4 );
    tmp( 2, 4 ) = y_off( i );
    anim.Robots{ i }.H_base = tmp;    
    robot.H_base = tmp;
    
    anim.attachRobot( robot )  
    
    N = N_moment( i );
    b = b_alpha( i );
    q = q_arr_demo( :, N );
    
    anim.Robots{ i }.updateKinematics( q );

    for j = 1 : 8
        if ismember( j , 1:8 )
            set( anim.gPatches{ i }{ j }, 'FaceAlpha', 0.2 )
        else
            set( anim.gPatches{ i }{ j }, 'FaceAlpha', b )
        end
    end
    
    y = y_off( i );
    scatter3( anim.hAxes,x+ p_arr_raw( 1, N), y + p_arr_raw( 2,N ), p_arr_raw( 3,N ), 100, 'filled', 'o', 'markerfacecolor', 'white', 'markeredgecolor', [0,0,0], 'linewidth', 5)
    
    xtmp = x+ p_arr_raw( 1, N );
    ytmp = y+ p_arr_raw( 2, N );
    ztmp =    p_arr_raw( 3, N );
    
    r1 = scl * R_arr_raw( :, 1, N );
    r2 = scl * R_arr_raw( :, 2, N );
    r3 = scl * R_arr_raw( :, 3, N );
    quiver3( anim.hAxes, xtmp, ytmp, ztmp, r1( 1 ), r1( 2 ), r1( 3 ), 'linewidth', 7, 'color', 'r' );
    quiver3( anim.hAxes, xtmp, ytmp, ztmp, r2( 1 ), r2( 2 ), r2( 3 ), 'linewidth', 7, 'color', 'g' );
    quiver3( anim.hAxes, xtmp, ytmp, ztmp, r3( 1 ), r3( 2 ), r3( 3 ), 'linewidth', 7, 'color', 'b' );
    
  
end

anim.update( 0 );
mySaveFig( gcf, 'revision_MIM_pos_and_orient_KUKA_iiwa14' )
