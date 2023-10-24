%% ==================================================================
%% [Title] Imitation Learning, in 3D Space, Position
% Author: Moses Chong-ook Nah
%  Email: mosesnah@mit.edu
%   Date: 2023.08.15
%% ==================================================================

%% [0A] Initialization
%% Initialization and call of data
clear all; close all; clc;
file_name = './data/example7/data1';
raw_data = parse_txt( [ file_name, '.txt' ], 0 );

t_arr = raw_data( :, 1 )' - raw_data( 1, 1 );
q_arr = raw_data( :, 2:8 )';
p0_arr = raw_data( :, 13:end )';

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

data1.y_arr = data1.y_arr + p0_arr( :, 1 );
data2.y_arr = data2.y_arr + data1.y_arr( :, end ) + [ 0.0,  0.1, -0.1 ]';
data3.y_arr = data3.y_arr + data2.y_arr( :, end ) + [ 0.0, -0.4,  0.1 ]';
data4.y_arr = data4.y_arr + data3.y_arr( :, end ) + [ 0.0,  0.3, -0.1 ]';

c_arr = [ 0.0000, 0.4470, 0.7410;
          0.8500, 0.3250, 0.0980;
          0.4940, 0.1840, 0.5560;
          0.4660, 0.6740, 0.1880 ];



%% ---- [1A] Call iiwa14 using Explicit-MATLAB

% Set figure size and attach robot to simulation
robot = iiwa14( 'high' );
robot.init( );

anim = Animation( 'Dimension', 3, 'xLim', [-0.1,1.1], 'yLim', [-0.7,0.5], 'zLim', [0,1.2], 'isSaveVideo', true );
anim.init( );
anim.attachRobot( robot )  
    
plot3(    anim.hAxes, p0_arr( 1, :   ), p0_arr( 2, :   ), p0_arr( 3, :   ), 'linewidth', 3, 'linestyle', '--', 'color', 'k' );
% plot3(    anim.hAxes, p_arr( 1, :   ), p_arr( 2, :   ), p_arr( 3, :   ), 'linewidth', 3, 'linestyle', '-', 'color', [0.8500, 0.3250, 0.0980]	);

view([90,0])

plot3( anim.hAxes, data1.y_arr( 1, : ), data1.y_arr( 2, : ), data1.y_arr( 3, : ), 'color', c_arr( 1, : ), 'linewidth', 5 )
plot3( anim.hAxes, data2.y_arr( 1, : ), data2.y_arr( 2, : ), data2.y_arr( 3, : ), 'color', c_arr( 2, : ), 'linewidth', 5 )
plot3( anim.hAxes, data3.y_arr( 1, : ), data3.y_arr( 2, : ), data3.y_arr( 3, : ), 'color', c_arr( 3, : ), 'linewidth', 5 )
plot3( anim.hAxes, data4.y_arr( 1, : ), data4.y_arr( 2, : ), data4.y_arr( 3, : ), 'color', c_arr( 4, : ), 'linewidth', 5 )

mk = 100;
scatter3( anim.hAxes, data1.y_arr( 1, 1 ), data1.y_arr( 2, 1 ), data1.y_arr( 3, 1 ), mk , 'filled', 'o'     , 'markerfacecolor', c_arr( 1, : ), 'markeredgecolor', 'black' );
scatter3( anim.hAxes, data1.y_arr( 1, end ), data1.y_arr( 2, end ), data1.y_arr( 3, end ), mk , 'filled', 'square', 'markerfacecolor', c_arr( 1, : ), 'markeredgecolor', 'black' );

scatter3( anim.hAxes, data2.y_arr( 1, 1 ), data2.y_arr( 2, 1 ), data2.y_arr( 3, 1 ), mk , 'filled', 'o'     , 'markerfacecolor', c_arr( 2, : ), 'markeredgecolor', 'black' );
scatter3( anim.hAxes, data2.y_arr( 1, end ), data2.y_arr( 2, end ), data2.y_arr( 3, end ), mk , 'filled', 'square', 'markerfacecolor', c_arr( 2, : ), 'markeredgecolor', 'black' );

scatter3( anim.hAxes, data3.y_arr( 1, 1 ), data3.y_arr( 2, 1 ), data3.y_arr( 3, 1 ), mk , 'filled', 'o'     , 'markerfacecolor', c_arr( 3, : ), 'markeredgecolor', 'black' );
scatter3( anim.hAxes, data3.y_arr( 1, end ), data3.y_arr( 2, end ), data3.y_arr( 3, end ), mk , 'filled', 'square', 'markerfacecolor', c_arr( 3, : ), 'markeredgecolor', 'black' );

scatter3( anim.hAxes, data4.y_arr( 1, 1 ), data4.y_arr( 2, 1 ), data4.y_arr( 3, 1 ), mk , 'filled', 'o'     , 'markerfacecolor', c_arr( 4, : ), 'markeredgecolor', 'black' );
scatter3( anim.hAxes, data4.y_arr( 1, end ), data4.y_arr( 2, end ), data4.y_arr( 3, end ), mk , 'filled', 'square', 'markerfacecolor', c_arr( 4, : ), 'markeredgecolor', 'black' );



for i = 1 : length( t_arr )
    robot.updateKinematics( q_arr( :, i ) );
    H = robot.getForwardKinematics( q_arr( :, i ) );
    anim.update( t_arr( i ) );
end

anim.close( );