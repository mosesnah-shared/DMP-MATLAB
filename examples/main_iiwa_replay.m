% [Title]     Replaying the iiwa from the joint trajectory
% [Author]    Moses Chong-ook Nah
% [Email]     mosesnah@mit.edu
% [Update]    At 2023.11.15

%% (--) Initialization
clear; close all; clc;
fig_config( 'fontSize', 20, 'markerSize', 10 )

%% =======================================================
%% (1-) Replaying the iiwa
%%  -- (1A) Call the Data

% Datasets are saved under data folder
file_name = './data/example6/test2';
raw_data  = parse_txt( [ file_name, '.txt' ], 0 );

% Time array, joint-trajectories, and p0 values
t_arr  = raw_data( :, 1 )' - raw_data( 1, 1 );
q_arr  = raw_data( :, 2:8 )';

% Need to comment out the p0 if the dataset do not include
% p0_arr = raw_data( :, 13:end )';


%%  -- (1B) Replay the Data

% Set figure size and attach robot to simulation
robot = iiwa14( 'high' );
robot.init( );

anim = Animation( 'Dimension', 3, 'xLim', [-0.1,1.1], 'yLim', [-0.7,0.5], 'zLim', [0,1.2], 'isSaveVideo', false );
anim.init( );
anim.attachRobot( robot )  
    
% p0 trajectory
% plot3(    anim.hAxes, p0_arr( 1, :   ), p0_arr( 2, :   ), p0_arr( 3, :   ), 'linewidth', 3, 'linestyle', '--', 'color', 'k' );

for i = 1 : length( t_arr )
    robot.updateKinematics( q_arr( :, i ) );
    H = robot.getForwardKinematics( q_arr( :, i ) );
    anim.update( t_arr( i ) );
end

anim.close( );