%% [Example Script] Data Filtering Script
% [Author] Moses Chong-ook Nah
% [Date]   2023.08.25

%% Initialization and call of data
clear all; close all; clc;
file_name = './data/example6/test3';
raw_data = parse_txt( [ file_name, '.txt' ], 0 );
t_arr = raw_data( :, 1 )' - raw_data( 1, 1 );
q_arr = raw_data( :, 2:8 )';

%% Call the robot without animation

% Set figure size and attach robot to simulation
robot = iiwa14( 'high' );
robot.init( );

%% Get the Forward Kinematics (p, R) data from text file.

[ ~, N ] = size( q_arr );

% End-effect Forward Kinematics, both position and orientation
p_arr_raw = zeros( 3, N );

% One should check whether R_arr is better with filtered q_arr
R_arr_raw = zeros( 3, 3, N );

for i = 1 : N
    robot.updateKinematics( q_arr( :, i ) );
    H = robot.getForwardKinematics( q_arr( :, i ) );
    p_arr_raw( :, i )    = H( 1:3, 4 );
    R_arr_raw( :, :, i ) = H( 1:3, 1:3 ); 
end

% Saving also the quaternion data
quat_arr_raw = zeros( 4, N );

for i = 1 : N
    quat = SO3_to_quat( R_arr_raw( :, :, i ) ); 
    quat_arr_raw( :, i ) = quat';
end

%% Filter out the p_arr 

dp_arr_raw  = data_diff(  p_arr_raw );
ddp_arr_raw = data_diff( dp_arr_raw );

% Filtering the Data
  p_arr_filt = zeros( size( p_arr_raw ), 'like', p_arr_raw );
 dp_arr_filt = zeros( size( p_arr_raw ), 'like', p_arr_raw );
ddp_arr_filt = zeros( size( p_arr_raw ), 'like', p_arr_raw );

for i = 1 : 3
   p_arr_filt( i, : ) = smoothdata( p_arr_raw( i, : ), "gaussian", 50 );

   tmp1 = data_diff( p_arr_filt );
   dp_arr_filt( i, : ) = smoothdata( tmp1( i, : ), "gaussian", 50 );

   tmp2 = data_diff( dp_arr_filt );
   ddp_arr_filt( i, : ) = smoothdata( tmp2( i, : ), "gaussian", 50 );   
end

tmp_plot = cell( 1, 6 );
cell{ 1 } = p_arr_raw;
cell{ 2 } = p_arr_filt;

cell{ 3 } = dp_arr_raw;
cell{ 4 } = dp_arr_filt;

cell{ 5 } = ddp_arr_raw;
cell{ 6 } = ddp_arr_filt;

tmp_label = [ "Pos.", "Vel.", "Acc." ];
f = figure( );
for i = 1:3
    for j = 1 :2
        subplot( 3, 2, 2*(i-1)+j )
        hold on
        plot( t_arr, cell{ 2*(i-1)+j }, 'linewidth', 3 )
        set( gca, 'xlim', [0, max( t_arr ) ], 'xticklabel', {}, 'fontsize', 30 )
        
        if j == 1
           ylabel( tmp_label( i ) )
           tmp = get( gca, 'ylim' );
        else
           set( gca, 'ylim', tmp ) 
        end
    end
end


%% Filter out the w_arr, dw_arr 

for i = 1 : N
    tmp = robot.getHybridJacobian( q_arr( :, i ) );
    JH_arr( :, :, i )  = tmp( 4:6, : );
end

% The raw w_arr 
 w_arr_raw  = zeros( size( p_arr_raw ), 'like', p_arr_raw );
dw_arr_raw  = zeros( size( p_arr_raw ), 'like', p_arr_raw );
 w_arr_filt = zeros( size( p_arr_raw ), 'like', p_arr_raw );
dw_arr_filt = zeros( size( p_arr_raw ), 'like', p_arr_raw );

% The non-filtered dq_arr
dq_arr_raw  = data_diff( q_arr );
dq_arr_filt = zeros( size( q_arr ), 'like', q_arr );

for i = 1 : 7
    dq_arr_filt( i, : ) = smoothdata( dq_arr_raw( i, : ), "gaussian", 50 );
end

for i = 1 : N
     w_arr_raw( :, i ) = JH_arr( :, :, i ) * dq_arr_raw( :, i );
    w_arr_filt( :, i ) = JH_arr( :, :, i ) * dq_arr_filt( :, i );
end

 dw_arr_raw = data_diff( w_arr_raw  );
dw_arr_filt = data_diff( w_arr_filt );

for i = 1 : 3
    dw_arr_filt( i, : ) = smoothdata( dw_arr_filt( i, : ), "gaussian", 50 );
end
f = figure( );
subplot( 2, 2, 1 )
plot( t_arr, w_arr_raw, 'linewidth', 3 )
set( gca, 'xlim', [0, max( t_arr ) ] )
ylim1 = get( gca, 'ylim' );
ylabel( 'Ang. Vel. ($rad/s$)')

subplot( 2, 2, 2 )
plot( t_arr, w_arr_filt, 'linewidth', 3 )
set( gca, 'xlim', [0, max( t_arr ) ], 'ylim', ylim1 )

subplot( 2, 2, 3 )
plot( t_arr, dw_arr_raw, 'linewidth', 3 )
set( gca, 'xlim', [0, max( t_arr ) ], 'ytick', [-5, 0, 5, 10]*10^-5 )
ylim2 = get( gca, 'ylim' );
ylabel( 'Ang. Acc. ($rad/s^2$)')
xlabel( 'Time ($s$)')

subplot( 2, 2, 4 )
plot( t_arr, dw_arr_filt, 'linewidth', 3 )
set( gca, 'xlim', [0, max( t_arr ) ], 'ylim', ylim2, 'ytick', [-5, 0, 5, 10]*10^-5 )
xlabel( 'Time ($s$)')

%% Check for sign flip of the quaternion

tmp = diff( quat_arr_raw' )';

% If there exists a value that changes so drastically, find the indexes 
[r, c_arr, ] = find( abs( tmp ) >= 1 );

% If not empty
if( ~isempty( c_arr ) ) 
    % c is the index, flip the sign 
    for c = c_arr
        quat_arr_raw( :,c+1:end ) = -quat_arr_raw( :,c+1:end );
    end
end

%% Save/Export Data

% Saving the crucial data 
save( [file_name, '.mat'], 'w_arr_raw', 'w_arr_filt', 'dw_arr_raw',  'dw_arr_filt', ...
                           'p_arr_raw', 'p_arr_filt', 'dp_arr_raw', 'dp_arr_filt', 'ddp_arr_raw', 'ddp_arr_filt', ... 
                           'R_arr_raw', 'quat_arr_raw', 'q_arr', 't_arr' )
