%% ==================================================================
%% [Title] Imitation Learning, in 3D Space, Position
% Author: Moses Chong-ook Nah
%  Email: mosesnah@mit.edu
%   Date: 2023.08.15
%% ==================================================================

%% [0A] Initialization
%% Initialization and call of data
clear all; close all; clc;
file_name = './data/example4/draw_M_actual_data1';
raw_data = parse_txt( [ file_name, '.txt' ] );

t_arr = raw_data( :, 1 )' - raw_data( 1, 1 );
q_arr = raw_data( :, 2:8 )';
p_arr = raw_data( :, 9:11 )';
p0_arr = raw_data( :, 13:end )';



%% ---- [1A] Call iiwa14 using Explicit-MATLAB

% Set figure size and attach robot to simulation
robot = iiwa14( 'high' );
robot.init( );

anim = Animation( 'Dimension', 3, 'xLim', [-0.7,0.7], 'yLim', [-0.7,0.7], 'zLim', [0,1.4], 'isSaveVideo', true );
anim.init( );
anim.attachRobot( robot )  
    
plot3(    anim.hAxes, p0_arr( 1, :   ), p0_arr( 2, :   ), p0_arr( 3, :   ), 'linewidth', 3, 'linestyle', '--', 'color', 'k' );
plot3(    anim.hAxes, p_arr( 1, :   ), p_arr( 2, :   ), p_arr( 3, :   ), 'linewidth', 3, 'linestyle', '-', 'color', [0.8500, 0.3250, 0.0980]	);

view([90,0])


for i = 1 : length( t_arr )
    robot.updateKinematics( q_arr( :, i ) );
    H = robot.getForwardKinematics( q_arr( :, i ) );
    anim.update( t_arr( i ) );
end

anim.close( );