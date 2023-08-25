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

dp_arr_filt  = diff_w_filter(       p_arr, "gaussian", 50 );
ddp_arr_filt = diff_w_filter( dp_arr_filt, "gaussian", 50 );

dp_arr_raw   = diff_w_filter(       p_arr, "none", 50 );
ddp_arr_raw  = diff_w_filter(  dp_arr_raw, "none", 50 );


%%
% Get also the analytical form of the Jacobian MAtrix
q_arr_sym = sym('q',[1 7]);

JH_sym = robot.getHybridJacobian( q_arr_sym );

dJH_sym = 0;
for i = 1 : 7
    dJH_sym = dJH_sym + diff( JH_sym, q_arr_sym( i ) );
end

JH_sym  = simplify( JH_sym );
dJH_sym = simplify( JH_sym );

% Trim off the JH part
JH_sym  =  JH_sym( 4:6, : );
dJH_sym = dJH_sym( 4:6, : );

% Get the values 
JH_arr  = zeros( 3, 7, N );

for i = 1 : N
    JH_arr( :, :, i )  = double( subs( JH_sym, q_arr_sym, q_arr( :, i )' ) );
end


%% Filter the p_arr data to get p_arr_filt, dq_arr_filt, ddq_arr_filt

p_arr_filt   = zeros( size( p_arr ) ,'like', p_arr );
dp_arr_filt  = zeros( size( p_arr ) ,'like', p_arr );
ddp_arr_filt = zeros( size( p_arr ) ,'like', p_arr );

for i = 1 : 3
   p_arr_filt( i, : ) = smoothdata( p_arr( i, : ), "gaussian", 50 );
end

color_arr = [ 0.0000, 0.4470, 0.7410;  ...
              0.8500, 0.3250, 0.0980;  ...
              0.9290, 0.6940, 0.1250];	

f = figure( ); a = axes( 'parent', f );
for i = 1 : 3
    subplot( 3, 1, i )
    hold on
    plot( t_arr, p_arr( i, : ), 'linestyle', '-', 'color', color_arr( i, : ), 'linewidth', 3 )
%     plot( t_arr, p_arr_filt( i, : ), 'linestyle', '--', 'color', color_arr( i, : ), 'linewidth', 3 )
    set( gca, 'xlim', [0, max( t_arr ) ] )
end

% Calculate dp_arr_filt, ddp_arr_filt
dp_arr = diff( p_arr_filt, 1, 2 )./diff( t_arr );
dp_arr( :, end+1 ) = dp_arr( :, end);

% Filter this one more time
for i = 1 : 3
   dp_arr_filt( i, : ) = smoothdata( dp_arr( i, : ), "gaussian", 50 );
end

f = figure( ); a = axes( 'parent', f );
for i = 1 : 3
    subplot( 3, 1, i )
    hold on
    plot( t_arr, dp_arr( i, : ), 'linestyle','-', 'color', color_arr( i, : ), 'linewidth', 3 )
    plot( t_arr, dp_arr_filt( i, : ), 'linestyle', '--', 'color', color_arr( i, : ), 'linewidth', 3 )    
    set( gca, 'xlim', [0, max( t_arr ) ] )
end

% Calculate dp_arr_filt, ddp_arr_filt
ddp_arr = diff( dp_arr, 1, 2 )./diff( t_arr );
ddp_arr( :, end+1 ) = ddp_arr( :, end);

% Filter this one more time
for i = 1 : 3
   ddp_arr_filt( i, : ) = smoothdata( ddp_arr( i, : ), "gaussian", 50 );
end

f = figure( ); a = axes( 'parent', f );
for i = 1 : 3
    subplot( 3, 1, i )
    hold on
    plot( t_arr, ddp_arr( i, : ), 'linestyle','-', 'color', color_arr( i, : ), 'linewidth', 3 )
    plot( t_arr, ddp_arr_filt( i, : ), 'linestyle', '--', 'color', color_arr( i, : ), 'linewidth', 3 )        
    set( gca, 'xlim', [0, max( t_arr ) ] )
end



%% Calculate the w_arr, dw_arr arrays and filter

% We use the fact that w_arr  = Jr(q) dq_arr
%                      dq_arr = d Jr(q) dq_arr + Jr(q) ddq_arr

% First, filter q_arr, dq_arr and ddq_arr 
dq_arr  = zeros( size( q_arr ) ,'like', q_arr );
ddq_arr = zeros( size( q_arr ) ,'like', q_arr );

q_arr_filt   = zeros( size( q_arr ) ,'like', q_arr );
dq_arr_filt  = zeros( size( q_arr ) ,'like', q_arr );
ddq_arr_filt = zeros( size( q_arr ) ,'like', q_arr );

for i = 1 : 7
   q_arr_filt( i, : ) = smoothdata( q_arr( i, : ), "gaussian", 50 );
end

f = figure( ); a = axes( 'parent', f );
hold on
plot( t_arr, q_arr, 'linestyle', '-', 'linewidth', 3 )
plot( t_arr, q_arr_filt, 'linestyle', '--', 'linewidth', 3, 'color', 'k' )
set( gca, 'xlim', [0, max( t_arr ) ] )

% Calculate dp_arr_filt, ddp_arr_filt
dq_arr = diff( q_arr_filt, 1, 2 )./diff( t_arr );
dq_arr( :, end+1 ) = dq_arr( :, end);

% Filter this one more time
for i = 1 : 7
   dq_arr_filt( i, : ) = smoothdata( dq_arr( i, : ), "gaussian", 50 );
end

f = figure( ); a = axes( 'parent', f );
hold on
plot( t_arr, dq_arr, 'linestyle', '-', 'linewidth', 3 )
plot( t_arr, dq_arr_filt, 'linestyle', '--', 'linewidth', 3, 'color', 'k' ) 
    set( gca, 'xlim', [0, max( t_arr ) ] )

% Calculate dp_arr_filt, ddp_arr_filt
ddq_arr = diff( dq_arr, 1, 2 )./diff( t_arr );
ddq_arr( :, end+1 ) = ddq_arr( :, end);

% Filter this one more time
for i = 1 : 7
   ddq_arr_filt( i, : ) = smoothdata( ddq_arr( i, : ), "gaussian", 50 );
end

f = figure( ); a = axes( 'parent', f );
hold on
plot( t_arr, ddq_arr, 'linestyle', '-', 'linewidth', 3 )
plot( t_arr, ddq_arr_filt, 'linestyle', '--', 'linewidth', 3, 'color', 'k' )     
set( gca, 'xlim', [0, max( t_arr ) ] )

%% Get w_arr 

 w_arr = zeros( size( p_arr), 'like', p_arr );
dw_arr = zeros( size( p_arr), 'like', p_arr );
 w_arr_filt = zeros( size( p_arr), 'like', p_arr );
dw_arr_filt = zeros( size( p_arr), 'like', p_arr );

for i = 1  : N
    w_arr( :, i ) = JH_arr( :, :, i ) * dq_arr_filt( :, i );
end

for i = 1 : 3
   w_arr_filt( i, : ) = smoothdata( w_arr( i, : ), "gaussian", 50 );
end

f = figure( ); a = axes( 'parent', f );
hold on
plot( t_arr, w_arr, 'linestyle', '-', 'linewidth', 3 )
plot( t_arr, w_arr_filt, 'linestyle', '--', 'linewidth', 3, 'color', 'k' ) 
    set( gca, 'xlim', [0, max( t_arr ) ] )

% Calculate dp_arr_filt, ddp_arr_filt
dw_arr = diff( w_arr_filt, 1, 2 )./diff( t_arr );
dw_arr( :, end+1 ) = dw_arr( :, end);

% Filter this one more time
for i = 1 : 3
   dw_arr_filt( i, : ) = smoothdata( dw_arr( i, : ), "gaussian", 50 );
end

f = figure( ); a = axes( 'parent', f );
hold on
plot( t_arr, dw_arr, 'linestyle', '-', 'linewidth', 3 )
plot( t_arr, dw_arr_filt, 'linestyle', '--', 'linewidth', 3, 'color', 'k' ) 
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