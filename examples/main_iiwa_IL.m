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


%% Get the Forward Kinematics (p, R) data from text file.

raw_data = parse_txt( './data/iiwa_example1.txt' );
t_arr = raw_data( :, 1 )' - raw_data( 1, 1 );
q_arr = raw_data( :, 2:end )';

[ ~, N ] = size( q_arr );

% End-effect Forward Kinematics, both position and orientation
p_arr = zeros( 3, N );

% One should check whether R_arr is better with filtered q_arr
R_arr = zeros( 3, 3, N );

for i = 1 : N
    robot.updateKinematics( q_arr( :, i ) );
    H = robot.getForwardKinematics( q_arr( :, i ) );
    p_arr( :, i )    = H( 1:3, 4 );
    R_arr( :, :, i ) = H( 1:3, 1:3 ); 
end

% Saving also the quaternion data
quat_data = zeros( 4, N );

for i = 1 : N
    quat = SO3_to_quat( R_arr( :, :, i ) ); 
    quat_data( :, i ) = quat';
end

%% Get the position, velocity and acceleration data of the end-effector. 

dp_arr  = data_diff(  p_arr );
ddp_arr = data_diff( dp_arr );

subplot( 3, 2, 1 )
hold on
plot( t_arr, p_arr, 'linewidth', 3 )
set( gca, 'xlim', [0, max( t_arr ) ], 'xticklabel', {}, 'fontsize', 30 )
ylabel( 'Pos. ($m$)')
ylim1 = get( gca, 'ylim' );

subplot( 3, 2, 3 )
plot( t_arr, dp_arr, 'linewidth', 3 )
set( gca, 'xlim', [0, max( t_arr ) ], 'xticklabel', {}, 'fontsize', 30 )
ylabel( 'Vel. ($m/s$)')
ylim2 = get( gca, 'ylim' );

subplot( 3, 2, 5 )
plot( t_arr, ddp_arr, 'linewidth', 3 )
set( gca, 'xlim', [0, max( t_arr ) ], 'fontsize', 30 )
ylabel( 'Acc. ($m/s^2$)')
xlabel( 'Time ($s$)')
legend( 'X', 'Y', 'Z', 'fontsize', 25 )
ylim3 = get( gca, 'ylim' );

% Filtering the Data
  p_arr_filt = zeros( size( p_arr ), 'like', p_arr );
 dp_arr_filt = zeros( size( p_arr ), 'like', p_arr );
ddp_arr_filt = zeros( size( p_arr ), 'like', p_arr );

for i = 1 : 3
   p_arr_filt( i, : ) = smoothdata( p_arr( i, : ), "gaussian", 50 );

   tmp1 = data_diff( p_arr_filt );
   dp_arr_filt( i, : ) = smoothdata( tmp1( i, : ), "gaussian", 50 );

   tmp2 = data_diff( dp_arr_filt );
   ddp_arr_filt( i, : ) = smoothdata( tmp2( i, : ), "gaussian", 50 );   
end

subplot( 3, 2, 2 )
hold on
plot( t_arr, p_arr_filt, 'linewidth', 3 )
set( gca, 'xlim', [0, max( t_arr ) ], 'xticklabel', {}, 'fontsize', 30, 'ylim', ylim1 )
ylabel( 'Pos. ($m$)')

subplot( 3, 2, 4 )
plot( t_arr, dp_arr_filt, 'linewidth', 3 )
set( gca, 'xlim', [0, max( t_arr ) ], 'xticklabel', {}, 'fontsize', 30, 'ylim', ylim2 )
ylabel( 'Vel. ($m/s$)')

subplot( 3, 2, 6 )
plot( t_arr, ddp_arr_filt, 'linewidth', 3 )
set( gca, 'xlim', [0, max( t_arr ) ], 'fontsize', 30, 'ylim', ylim3 )
ylabel( 'Acc. ($m/s^2$)')
xlabel( 'Time ($s$)')
legend( 'X', 'Y', 'Z', 'fontsize', 25 )



%% Get also the analytical form of the Jacobian MAtrix

for i = 1 : N
    tmp = robot.getHybridJacobian( q_arr( :, i ) );
    JH_arr( :, :, i )  = tmp( 4:6, : );
end

% The raw w_arr 
 w_arr = zeros( size( p_arr ), 'like', p_arr );
dw_arr = zeros( size( p_arr ), 'like', p_arr );

% The filtered w_arr 
 w_arr_filt = zeros( size( p_arr ), 'like', p_arr );
dw_arr_filt = zeros( size( p_arr ), 'like', p_arr );

% The non-filtered dq_arr
dq_arr = data_diff( q_arr );
dq_arr_filt = zeros( size( q_arr ), 'like', q_arr );

for i = 1 : 7
    dq_arr_filt( i, : ) = smoothdata( dq_arr( i, : ), "gaussian", 50 );
end

for i = 1 : N
         w_arr( :, i ) = JH_arr( :, :, i ) * dq_arr( :, i );
    w_arr_filt( :, i ) = JH_arr( :, :, i ) * dq_arr_filt( :, i );
end

 dw_arr     = data_diff( w_arr      );
dw_arr_filt = data_diff( w_arr_filt );

for i = 1 : 3
    dw_arr_filt( i, : ) = smoothdata( dw_arr_filt( i, : ), "gaussian", 50 );
end

subplot( 2, 2, 1 )
plot( t_arr, w_arr, 'linewidth', 3 )
set( gca, 'xlim', [0, max( t_arr ) ] )
ylim1 = get( gca, 'ylim' );

subplot( 2, 2, 2 )
plot( t_arr, w_arr_filt, 'linewidth', 3 )
set( gca, 'xlim', [0, max( t_arr ) ], 'ylim', ylim1 )

subplot( 2, 2, 3 )
plot( t_arr, dw_arr, 'linewidth', 3 )
set( gca, 'xlim', [0, max( t_arr ) ] )
ylim2 = get( gca, 'ylim' )

subplot( 2, 2, 4 )
plot( t_arr, dw_arr_filt, 'linewidth', 3 )
set( gca, 'xlim', [0, max( t_arr ) ], 'ylim', ylim2 )


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