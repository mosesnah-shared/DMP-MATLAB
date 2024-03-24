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
% file_name = './data/example4/iiwa_example_pos';
% file_name = './data/example3/iiwa_example_orient';
file_name = './data/example8/modular_input';
raw_data  = parse_txt( [ file_name, '.txt' ], 0 );

% Time array, joint-trajectories, and p0 values
t_arr  = raw_data( :, 1 )' - raw_data( 1, 1 );
q_arr  = raw_data( :, 2:8 )';

% Set figure size and attach robot to simulation
robot = iiwa14( 'high' );
robot.init( );

Nt = length( t_arr );
p_arr = zeros( 3, Nt );
R_arr = zeros( 3, 3, Nt );
for i = 1 : Nt
    q = q_arr( :, i );

    H = robot.getForwardKinematics( q );
    p_arr( :, i ) = H(1:3,4 ); 
    R_arr( :, :, i ) = H(1:3,1:3 ); 
end


%%  -- (1B) Replay the Data

anim = Animation( 'Dimension', 3, 'xLim', [-0.1,1.1], 'yLim', [-0.6,0.6], 'zLim', [0,1.2], 'isSaveVideo', true );
anim.init( );
anim.attachRobot( robot )  
    
% p0 trajectory
poff = 0.2;
scatter3( anim.hAxes, p_arr( 1, 1 )+2*poff , p_arr( 2, 1 ), p_arr( 3, 1 ), 400, 'filled', 'o', 'markerfacecolor', 'w', 'markeredgecolor', [0.6350 0.0780 0.1840], 'linewidth', 5 );
scatter3( anim.hAxes, p_arr( 1, end )+2*poff , p_arr( 2, end ), p_arr( 3, end ), 400, 'filled', 'd', 'markerfacecolor', 'w', 'markeredgecolor', [0.6350 0.0780 0.1840], 'linewidth', 5 );

scl = 0.1;
R1 = scl*R_arr( :, :, 1   );
R2 = scl*R_arr( :, :, end );

quiver3( anim.hAxes, p_arr( 1, 1 )+2*poff , p_arr( 2, 1 ), p_arr( 3, 1 ), R1( 1,1 ),R1( 2,1 ), R1( 3,1 ), 'linewidth', 5, 'color', 'r')
quiver3( anim.hAxes, p_arr( 1, 1 )+2*poff , p_arr( 2, 1 ), p_arr( 3, 1 ), R1( 1,2 ),R1( 2,2 ), R1( 3,2 ), 'linewidth', 5, 'color', 'g')
quiver3( anim.hAxes, p_arr( 1, 1 )+2*poff , p_arr( 2, 1 ), p_arr( 3, 1 ), R1( 1,3 ),R1( 2,3 ), R1( 3,3 ), 'linewidth', 5, 'color', 'b')

gtmp1 = quiver3( anim.hAxes, p_arr( 1, 1 )+2*poff , p_arr( 2, 1 ), p_arr( 3, 1 ), R1( 1,1 ),R1( 2,1 ), R1( 3,1 ), 'linewidth', 5, 'color', 'r');
gtmp2 = quiver3( anim.hAxes, p_arr( 1, 1 )+2*poff , p_arr( 2, 1 ), p_arr( 3, 1 ), R1( 1,2 ),R1( 2,2 ), R1( 3,2 ), 'linewidth', 5, 'color', 'g');
gtmp3 = quiver3( anim.hAxes, p_arr( 1, 1 )+2*poff , p_arr( 2, 1 ), p_arr( 3, 1 ), R1( 1,3 ),R1( 2,3 ), R1( 3,3 ), 'linewidth', 5, 'color', 'b');

quiver3( anim.hAxes, p_arr( 1, end )+2*poff , p_arr( 2, end ), p_arr( 3, end ), R1( 1,1 ),R1( 2,1 ), R1( 3,1 ), 'linewidth', 5, 'color', 'r')
quiver3( anim.hAxes, p_arr( 1, end )+2*poff , p_arr( 2, end ), p_arr( 3, end ), R1( 1,2 ),R1( 2,2 ), R1( 3,2 ), 'linewidth', 5, 'color', 'g')
quiver3( anim.hAxes, p_arr( 1, end )+2*poff , p_arr( 2, end ), p_arr( 3, end ), R1( 1,3 ),R1( 2,3 ), R1( 3,3 ), 'linewidth', 5, 'color', 'b')

N1 = 2800;
R1 = scl*R_arr( :, :, N1 );
scatter3( anim.hAxes, p_arr( 1, N1 )+2*poff , p_arr( 2, N1 ), p_arr( 3, N1 ), 400, 'filled', 'o', 'markerfacecolor', 'w', 'markeredgecolor', [0.6350 0.0780 0.1840], 'linewidth', 5 );
quiver3( anim.hAxes, p_arr( 1, N1 )+2*poff , p_arr( 2, N1 ), p_arr( 3, N1 ), R1( 1,1 ),R1( 2,1 ), R1( 3,1 ), 'linewidth', 5, 'color', 'r')
quiver3( anim.hAxes, p_arr( 1, N1 )+2*poff , p_arr( 2, N1 ), p_arr( 3, N1 ), R1( 1,2 ),R1( 2,2 ), R1( 3,2 ), 'linewidth', 5, 'color', 'g')
quiver3( anim.hAxes, p_arr( 1, N1 )+2*poff , p_arr( 2, N1 ), p_arr( 3, N1 ), R1( 1,3 ),R1( 2,3 ), R1( 3,3 ), 'linewidth', 5, 'color', 'b')

N1 = 3300;
R1 = scl*R_arr( :, :, N1 );
scatter3( anim.hAxes, p_arr( 1, N1 )+2*poff , p_arr( 2, N1 ), p_arr( 3, N1 ), 400, 'filled', 'o', 'markerfacecolor', 'w', 'markeredgecolor', [0.6350 0.0780 0.1840], 'linewidth', 5 );
quiver3( anim.hAxes, p_arr( 1, N1 )+2*poff , p_arr( 2, N1 ), p_arr( 3, N1 ), R1( 1,1 ),R1( 2,1 ), R1( 3,1 ), 'linewidth', 5, 'color', 'r')
quiver3( anim.hAxes, p_arr( 1, N1 )+2*poff , p_arr( 2, N1 ), p_arr( 3, N1 ), R1( 1,2 ),R1( 2,2 ), R1( 3,2 ), 'linewidth', 5, 'color', 'g')
quiver3( anim.hAxes, p_arr( 1, N1 )+2*poff , p_arr( 2, N1 ), p_arr( 3, N1 ), R1( 1,3 ),R1( 2,3 ), R1( 3,3 ), 'linewidth', 5, 'color', 'b')

N1 = 3800;
R1 = scl*R_arr( :, :, N1 );
scatter3( anim.hAxes, p_arr( 1, N1 )+2*poff , p_arr( 2, N1 ), p_arr( 3, N1 ), 400, 'filled', 'o', 'markerfacecolor', 'w', 'markeredgecolor', [0.6350 0.0780 0.1840], 'linewidth', 5 );
quiver3( anim.hAxes, p_arr( 1, N1 )+2*poff , p_arr( 2, N1 ), p_arr( 3, N1 ), R1( 1,1 ),R1( 2,1 ), R1( 3,1 ), 'linewidth', 5, 'color', 'r')
quiver3( anim.hAxes, p_arr( 1, N1 )+2*poff , p_arr( 2, N1 ), p_arr( 3, N1 ), R1( 1,2 ),R1( 2,2 ), R1( 3,2 ), 'linewidth', 5, 'color', 'g')
quiver3( anim.hAxes, p_arr( 1, N1 )+2*poff , p_arr( 2, N1 ), p_arr( 3, N1 ), R1( 1,3 ),R1( 2,3 ), R1( 3,3 ), 'linewidth', 5, 'color', 'b')

plot3(    anim.hAxes, p_arr( 1, : )+poff, p_arr( 2, : ), p_arr( 3, :   ), 'linewidth', 6, 'linestyle', '-', 'color',[0.6350 0.0780 0.1840]);

ptmp = scatter3( anim.hAxes, p_arr( 1, 1 )+2*poff , p_arr( 2, 1 ), p_arr( 3, 1 ), 200, 'filled', 'o', 'markerfacecolor', 'w', 'markeredgecolor', [0.6350 0.0780 0.1840], 'linewidth', 5 );

set( anim.hAxes, 'xticklabel', {}, 'yticklabel', {}, 'zticklabel', {} )
set( anim.SubTitle, 'visible', 'off' )
set( anim.SubTitle, 'visible', 'off' )
set( anim.hAxes.XLabel, 'visible', 'off' )
set( anim.hAxes.YLabel, 'visible', 'off' )
set( anim.hAxes.ZLabel, 'visible', 'off' )

view( anim.hAxes, [ 90, 0 ])
for i = 1 : length( t_arr )
    robot.updateKinematics( q_arr( :, i ) );
    H = robot.getForwardKinematics( q_arr( :, i ) );
    set( ptmp, 'Xdata', p_arr( 1, i )+2*poff, 'Ydata', p_arr( 2, i ), 'Zdata', p_arr( 3, i ) )

    Rtmp = scl*R_arr( :, :,i );

    set( gtmp1, 'Xdata', p_arr( 1, i )+2*poff, 'Ydata', p_arr( 2, i ), 'Zdata', p_arr( 3, i ), 'Udata', Rtmp( 1, 1 ), 'VData', Rtmp( 2, 1 ), 'Wdata', Rtmp( 3, 1 )  )
    set( gtmp2, 'Xdata', p_arr( 1, i )+2*poff, 'Ydata', p_arr( 2, i ), 'Zdata', p_arr( 3, i ), 'Udata', Rtmp( 1, 2 ), 'VData', Rtmp( 2, 2 ), 'Wdata', Rtmp( 3, 2 )  )
    set( gtmp3, 'Xdata', p_arr( 1, i )+2*poff, 'Ydata', p_arr( 2, i ), 'Zdata', p_arr( 3, i ), 'Udata', Rtmp( 1, 3 ), 'VData', Rtmp( 2, 3 ), 'Wdata', Rtmp( 3, 3 )  )

    anim.update( t_arr( i ) );
end


anim.close( );