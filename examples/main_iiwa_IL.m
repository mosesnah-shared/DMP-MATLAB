%% [Example Script] Plotting the Nonlinear Forcing Terms
% [Author] Moses Chong-ook Nah
% [Date]   2023.08.15

%% [Initialization] 
clear; close all; clc;
fig_config( 'fontSize', 20, 'markerSize', 10 )

%% Call iiwa from Explicit-MATLAB

% Set figure size and attach robot to simulation
robot = iiwa14( 'high' );
robot.init( );

anim = Animation( 'Dimension', 3, 'xLim', [-0.7,0.7], 'yLim', [-0.7,0.7], 'zLim', [0,1.4], 'isSaveVideo', true );
anim.init( );
anim.attachRobot( robot )  

view( 90, 0 );


%% Get the data from txt

raw_data = parse_txt( './data/iiwa_example1.txt' );
t_arr = raw_data( :, 1 )' - raw_data( 1, 1 );
q_arr = raw_data( :, 2:end );

[ N, ~ ] = size( q_arr );

% End-effect Forward Kinematics, both position and orinetation
p_arr = zeros( 3, N );
R_arr = zeros( 3, 3, N );

for i = 1 : N
    robot.updateKinematics( q_arr( i, : ) );
    H = robot.getForwardKinematics( q_arr( i, : ) );
    p_arr( :, i )    = H( 1:3, 4 );
    R_arr( :, :, i ) = H( 1:3, 1:3 ); 
end

% Saving also the quaternion data
quat_data = zeros( 4, N );

for i = 1 : N
    quat = SO3_to_quat( R_arr( :, :, i ) ); 
    quat_data( :, i ) = quat';
end


%% Cleaning up the q_arr data

% Before Filtering 
subplot( 2,1,1 )
plot( t_arr, q_arr )
set( gca, 'xlim', [0, max( t_arr ) ] )

q_arr_filt = zeros( size( q_arr ) ,'like', q_arr );

for i = 1 : 7
   q_arr_filt( :, i ) = smoothdata( q_arr( :, i ), "gaussian" );
end

subplot( 2,1,2 )
plot( t_arr, q_arr_filt )
set( gca, 'xlim', [0, max( t_arr ) ] )

%% Filter the p_arr data to get p_arr_filt, dq_arr_filt, ddq_arr_filt
close all;
p_arr_filt = zeros( size( p_arr ) ,'like', p_arr );

for i = 1 : 3
   p_arr_filt( i, : ) = smoothdata( p_arr( i, : ), "gaussian", 100 );
end

color_arr = [ 0.0000, 0.4470, 0.7410;  ...
              0.8500, 0.3250, 0.0980;  ...
              0.9290, 0.6940, 0.1250];	

for i = 1 : 3
    f = figure( ); a = axes( 'parent' ,f );
    hold on
    plot( t_arr, p_arr( i, : ), 'linestyle', '-', 'color', color_arr( i, : ), 'linewidth', 3 )
    plot( t_arr, p_arr_filt( i, : ), 'linestyle', '--', 'color', color_arr( i, : ), 'linewidth', 3 )
    set( gca, 'xlim', [0, max( t_arr ) ] )
end


%% Calculate dq_arr, ddq_arr from data

dq_arr_filt = diff( q_arr_filt, 1, 1 );
dq_arr_filt( end + 1, : ) = dq_arr_filt( end, : );
subplot( 2, 1, 1 )
plot( t_arr, dq_arr_filt )
set( gca, 'xlim', [0, max( t_arr ) ] )

% Filtering again 

ddq_arr_filt = diff( dq_arr_filt, 1, 1 );
ddq_arr_filt( end + 1, : ) = ddq_arr_filt( end, : );
subplot( 2, 1, 2 )
plot( t_arr, ddq_arr_filt )
set( gca, 'xlim', [0, max( t_arr ) ] )

%% Running the Animation 

plot3( anim.hAxes, p_arr( 1, : ), p_arr( 2, : ), p_arr( 3, : ), 'linewidth', 3, 'linestyle', '--', 'color', 'k' );
scatter3( anim.hAxes, p_arr( 1, 1   ), p_arr( 2, 1   ), p_arr( 3, 1   ), 300, 'filled', 'linewidth', 3, 'markerfacecolor', [0 0.4470 0.7410], 'markeredgecolor', 'k'  );
scatter3( anim.hAxes, p_arr( 1, end ), p_arr( 2, end ), p_arr( 3, end ), 300, 'filled', 'linewidth', 3, 'markerfacecolor', [0 0.4470 0.7410], 'markeredgecolor', 'k' );

view([109.1242, 2.7291])

for i = 1 : N
    robot.updateKinematics( q_arr( i, : ) );
    H = robot.getForwardKinematics( q_arr( i, : ) );
    anim.update( t_arr( i ) );
end

anim.close( );