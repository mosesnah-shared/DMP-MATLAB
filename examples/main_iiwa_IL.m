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

% End-effect Kinematics
  p_arr = zeros( 3, N );
 dp_arr = zeros( 3, N );
ddp_arr = zeros( 3, N );

for i = 1 : N
    robot.updateKinematics( q_arr( i, : ) );
    H = robot.getForwardKinematics( q_arr( i, : ) );
    p_arr( :, i ) = H( 1:3, 4 );
    
end

%% Running the Animation 

plot3( anim.hAxes, p_arr( 1, : ), p_arr( 2, : ), p_arr( 3, : ), 'linewidth', 3, 'linestyle', '--', 'color', 'k' );
scatter3( anim.hAxes, p_arr( 1, 1   ), p_arr( 2, 1   ), p_arr( 3, 1   ), 300, 'filled', 'linewidth', 3, 'markerfacecolor', [0 0.4470 0.7410], 'markeredgecolor', 'k'  );
scatter3( anim.hAxes, p_arr( 1, end ), p_arr( 2, end ), p_arr( 3, end ), 300, 'filled', 'linewidth', 3, 'markerfacecolor', [0 0.4470 0.7410], 'markeredgecolor', 'k' );

for i = 1 : N
    robot.updateKinematics( q_arr( i, : ) );
    H = robot.getForwardKinematics( q_arr( i, : ) );
    anim.update( t_arr( i ) );
end



anim.close( );