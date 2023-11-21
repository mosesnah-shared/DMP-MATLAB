% [Title]     Main Contraction Theory Position and Orientation Example
% [Author]    Moses Chong-ook Nah
% [Email]     mosesnah@mit.edu
% [Update]    At 2023.11.16

%% (--) Initialization
clear; close all; clc;
fig_config( 'fontSize', 20, 'markerSize', 10 )

%% ===========================================================
%% (1-) Data Collection and Trajectory Generation for Position
%%  -- (1A) Calling the trajectory for position and orientation

% The string list for the types of trajectories
% We have two trajectories, both position and orientation
traj_names = [ "drawM" , "lift_up_down" ];

Ntraj = length( traj_names );

traj_pos_data    = cell( 1, Ntraj );
traj_orient_data = cell( 1, Ntraj );

% Loading the data
for i = 1 : Ntraj
    % For position
    tmp = load( [ './learned_parameters/', traj_names{ i }, '_pos.mat' ] );
    traj_pos_data{ i } = tmp.data;

    % For orientation
    tmp = load( [ './learned_parameters/', traj_names{ i }, '_orient.mat' ] );
    traj_orient_data{ i } = tmp.data;
end

% Set as same duration
traj_pos_data{ 1 }.tau    = 3;
traj_pos_data{ 2 }.tau    = 3;
traj_orient_data{ 1 }.tau = 3;
traj_orient_data{ 2 }.tau = 3;

% Color array for plots
c_arr = [ 0.0000, 0.4470, 0.7410;
          0.8500, 0.3250, 0.0980;
          0.9290, 0.6940, 0.1250;
          0.4940, 0.1840, 0.5560;
          0.4660, 0.6740, 0.1880;
          0.3010, 0.7450, 0.9330];

clear tmp*

%%  -- (1B) Generating the trajectories for Position 

tinit =    1.0;           % The initial time of the simulation
T     =    8.0;           % The   whole time of the simulation 
toff  =   -0.5;           % Time-offset for the subsequent movement to start, must be negative
dt    =   1e-3;           % Time-step  for the simulation
t_arr = 0:dt:T;           % Time array for the simulation

Nt = length( t_arr );

% Time offset must be strictly negative. 
assert( toff < 0 );

% Position and velocity array of the trajectory
 p_data_arr = zeros( 3, length( t_arr ), Ntraj );
dp_data_arr = zeros( 3, length( t_arr ), Ntraj );

% error vector and its velocity for orientation
 e_data_arr = zeros( 3, length( t_arr ), Ntraj );
de_data_arr = zeros( 3, length( t_arr ), Ntraj );

% Scaling and Rotation of the trajectory, in position
scl_arr = [ 2.5, 5.0 ];

R_arr = zeros( 3, 3, Ntraj );
R_arr( :, :, 1 ) =  rotx( 0 ) * roty( 0 ) * rotz( 0 );
R_arr( :, :, 2 ) =  rotx( 0 ) * roty( 0 ) * rotz( 0 );

% The input force array for the movements
pos_force_arr = zeros( 3, length( t_arr )-1, Ntraj );
pos_goal_arr = zeros( 3, Ntraj );

% The goal and initated time and final time for each movement
t0i_arr = zeros( 1, Ntraj );
t0f_arr = zeros( 1, Ntraj );

% The first movement starts at tinit, hence setting up
t0i = tinit;          

% We need to first generate the DMP 
for i = 1 : Ntraj
    
    % ================================================================== %
    % ========================= Position Data ========================== %
    % ================================================================== %
    % Call the data from the position trajectory
    data = traj_pos_data{ i };

    % Saving the initial time in the array
    t0i_arr( i ) = t0i;
    t0f_arr( i ) = t0i + data.tau;
    
    % Getting the number of basis functions from the nonlinear forcing term
    [ ~, N ] = size( data.weight );

    % The Three elements of DMP, position
    cs        = CanonicalSystem( 'discrete', data.tau, data.alpha_s );
    fs        = NonlinearForcingTerm( cs, N );
    trans_sys = TransformationSystem( data.alpha_z, data.beta_z, cs );

    % Calculate the nonlinear forcing term 
    % This can be diminishing
    input_arr = fs.calc_forcing_term( t_arr( 1:end-1 ), data.weight, t0i, eye( 3 ) );

    % Rotation and Scaling
     R  = R_arr( :, :, i );
    scl = scl_arr( i );
        
    % Calculate the new goal location and save
    new_goal = scl*R*data.goal;
    pos_goal_arr( :, i ) = new_goal;

    % Calculate the new force and save
    new_force = scl*R*input_arr; 
    pos_force_arr( :, :, i ) = new_force;

    % Rollout, note that initial condition is set to be zeros. 
    [ y_arr, ~, dy_arr ] = trans_sys.rollout( zeros( 3, 1 ), zeros( 3, 1 ), new_goal, new_force, t0i, t_arr  );    
    
     p_data_arr( :, :, i ) =  y_arr;
    dp_data_arr( :, :, i ) = dy_arr;

    % ================================================================== %
    % ======================= Orientation Data ========================= %
    % ================================================================== %
    % Call the data from the orientation trajectory
    data = traj_orient_data{ i };

    % Getting the number of basis functions from the nonlinear forcing term
    [ ~, N ] = size( data.weight );

    % The Three elements of DMP, quaternion
    cs        = CanonicalSystem( 'discrete', data.tau, data.alpha_s );
    trans_sys = TransformationSystemQuat( data.alpha_z, data.beta_z, cs );
    fs        = NonlinearForcingTerm( cs, N );
    
    % Scaling Matrix
    eq0 = get_quat_error( data.init, data.goal );
    Sr  = diag( eq0 );

    input_arr = fs.calc_forcing_term( t_arr( 1:end-1 ), data.weight, t0i, Sr  );
    
    [ eq_arr, z_arr, deq_arr ] = trans_sys.rollout( eq0, zeros( 3, 1 ), 1*input_arr, t0i, t_arr  );

     e_data_arr( :, :, i ) =  eq_arr;
    de_data_arr( :, :, i ) = deq_arr;

    % Define the next movement's starting time.
    t0i = t0i + toff + data.tau;

end

% Add position offset 
poff = zeros( 3, 3 );
poff( :, 1 ) = [  0.1; 0.5;  0.0 ];
p_data_arr( :, :, 2 ) = p_data_arr( :, :, 2 ) + pos_goal_arr( :, 1 ) + poff( :, 1 );

pos_goal_arr( :, 1 ) = p_data_arr( :, end, 1 );
pos_goal_arr( :, 2 ) = p_data_arr( :, end, 2 );

% Plotting the positions to double cehck
f = figure( ); a = axes( 'parent', f );
axis equal; view( 3 );
hold( a, 'on' )

for i = 1 : Ntraj
    plot3( a, p_data_arr( 1, :, i ), p_data_arr( 2, :, i ), p_data_arr( 3, :, i ), 'linewidth', 5, 'color', c_arr( i, : ) )
    scatter3( a, p_data_arr( 1,   1, i ), p_data_arr( 2,   1, i ), p_data_arr( 3,   1, i ), 'filled', 'o', 'markerfacecolor', 'white', 'markeredgecolor', c_arr( i, : ), 'linewidth', 4 )
    scatter3( a, p_data_arr( 1, end, i ), p_data_arr( 2, end, i ), p_data_arr( 3, end, i ), 'filled', 'd', 'markerfacecolor', 'white', 'markeredgecolor', c_arr( i, : ), 'linewidth', 4 )    
end

% Also Plotting the orientation 
% First, Reviving the orientations from quat to R arr
R_data_arr = zeros( 3, 3, Nt, Ntraj );

for i = 1 : Ntraj
    % First, extend the 3D array to quaternion
    eq_arr = e_data_arr( :, :, i );

    for j = 1 : Nt
        tmp1 = R3_to_quat( eq_arr( :,j ) );
        tmp2 = quat_mul( traj_orient_data{ i }.goal, quat_conj( ExpQuat( 0.5 * tmp1 ) )' );

        R_data_arr( :, :, j, i ) = quat_to_SO3( tmp2 );
    end
end

% Plotting the orientation
tmp_scl = 0.1;
tmpN    = 30;
tmparr  = round( linspace( 1, Nt, tmpN ) );

for i = 1 : Ntraj

    R_arr = R_data_arr( :, :, :, i );
    p_arr = p_data_arr( :, :, i );

    for j = tmparr
        scatter3( a, p_arr( 1, j ), p_arr( 2, j ), p_arr( 3, j ), 200, 'filled', 'o', 'markeredgecolor', 'black', 'markerfacecolor', 'white', 'linewidth', 5 )    
         quiver3( a, p_arr( 1, j ), p_arr( 2, j ), p_arr( 3, j ), tmp_scl * R_arr( 1, 1, j ), tmp_scl * R_arr( 2, 1, j ), tmp_scl * R_arr( 3, 1, j ), 'linewidth', 4, 'color', 'r' )
         quiver3( a, p_arr( 1, j ), p_arr( 2, j ), p_arr( 3, j ), tmp_scl * R_arr( 1, 2, j ), tmp_scl * R_arr( 2, 2, j ), tmp_scl * R_arr( 3, 2, j ), 'linewidth', 4, 'color', 'g' )
         quiver3( a, p_arr( 1, j ), p_arr( 2, j ), p_arr( 3, j ), tmp_scl * R_arr( 1, 3, j ), tmp_scl * R_arr( 2, 3, j ), tmp_scl * R_arr( 3, 3, j ), 'linewidth', 4, 'color', 'b' )
    end
end

%%  -- (1C) Use Contraction Theory

% Setting the Departure and Arrival Time 
t_depart_arr = t0f_arr( 1 ) + toff + 0.1;
t_arrive_arr = t0i_arr( Ntraj ) + 0.5;

% Define the A matrix and az, bz values
% All should be identical, 
% Tau can be arbitrarily chosen based on the proof
% Choosing one arbitrarily
tau = traj_pos_data{ 1 }.tau;
as  = traj_pos_data{ 1 }.alpha_s;
az  = traj_pos_data{ 1 }.alpha_z;
bz  = traj_pos_data{ 1 }.beta_z;

% A Matrix
Amat = 1/tau * [            -as,       zeros( 1, 3 ),  zeros( 1, 3 ); 
                  zeros( 3, 1 ),       zeros( 3, 3 ),       eye( 3 );
                  zeros( 3, 1 ), -az * bz * eye( 3 ), -az * eye( 3 )];

Cg = 1/tau * [ zeros( 4, 1 ); az * bz * pos_goal_arr( :, 1 )];

% The initial condition is identical to the first trajectory x0 value.
x0 = [ 1; p_data_arr( :, 1, 1 ); dp_data_arr( :, 1, 1 )*traj_pos_data{1}.tau];

% The coupled x array 
xc_arr = zeros( 7, length( t_arr ) );
xc_arr( :, 1 ) = x0;

xcurr = x0;
% Iterating through the x_arr

% Integrating over time. 
for i = 1 : length( t_arr )-1
    
    t = t_arr( i );

    if t <= tinit
        dx = zeros( 7, 1 );
    else
        % Add the force for the first movement
        dx  = Amat * xcurr + 1/traj_pos_data{ 1 }.tau * [ zeros( 4, 1 ); pos_force_arr( :, i, 1 ) ] + Cg;

        % Iterating through the movement with activation
        for j = 2 : Ntraj
            fp = pos_force_arr( :, i, j-1 );  % Force for previous movement
            fn = pos_force_arr( :, i, j   );  % Force for the next movement
            
            Dp = traj_pos_data{ j-1 }.tau;
            Dn = traj_pos_data{ j   }.tau;

            td   = t_depart_arr( j-1 );
            ta   = t_arrive_arr( j-1 );
        
            gp = 1/Dp * [ zeros( 4, 1 ); az * bz * pos_goal_arr( :, j-1 ) ];
            gn = 1/Dn * [ zeros( 4, 1 ); az * bz * pos_goal_arr( :, j   ) ];
            
            gain = clip_func( t, td, ta );
            dx  = dx +  gain * ( - 1/Dp * [ zeros( 4, 1 ); fp ] + 1/Dn * [ zeros( 4, 1 ); fn ] ) ...
                     +  gain * ( -gp + gn ); 
        end

    end

    xc_arr( :, i+1 ) = xc_arr( :, i ) + dx * dt;     
    xcurr = xc_arr( :, i+1 );

end


%%  -- (1D) Generating the Video    

% v = VideoWriter( 'video.mp4','MPEG-4' );
% v.FrameRate = 30;

% open( v );f
Nstep = round( 1/dt / 30 );

t1 = title( sprintf( 'Time %.3f s', 0 ) );

axis equal; view( 3 );
hold( a, 'on' )
s_arr = cell( 1, Ntraj );
   qx = cell( 1, Ntraj );
   qy = cell( 1, Ntraj );
   qz = cell( 1, Ntraj );


% Draw the position marker and the quiver plot
for i = 1: Ntraj
    s_arr{ i } = scatter3( a, p_data_arr( 1, 1, i ), p_data_arr( 2, 1, i ), p_data_arr( 3, 1, i ), 100, 'filled', 'o', 'markerfacecolor', c_arr( i, : ), 'markeredgecolor', c_arr( i, : ), 'linewidth', 4 );

    plot3( a, p_data_arr( 1, :, i ), p_data_arr( 2, :, i ), p_data_arr( 3, :, i ), 'linewidth', 5, 'color', c_arr( i, : ) )
    scatter3( a, p_data_arr( 1,   1, i ), p_data_arr( 2,   1, i ), p_data_arr( 3,   1, i ), 20, 'filled', 'o', 'markerfacecolor', 'white', 'markeredgecolor', c_arr( i, : ), 'linewidth', 4 )
    scatter3( a, p_data_arr( 1, end, i ), p_data_arr( 2, end, i ), p_data_arr( 3, end, i ), 20, 'filled', 'd', 'markerfacecolor', 'white', 'markeredgecolor', c_arr( i, : ), 'linewidth', 4 )  

    % The xyz frame
    qx{ i } = quiver3( a, p_data_arr( 1, 1, i ), p_data_arr( 2, 1, i ), p_data_arr( 3, 1, i ), 2*tmp_scl * R_data_arr( 1, 1, 1, i ), 2*tmp_scl * R_data_arr( 2, 1, 1, i ), 2*tmp_scl * R_data_arr( 3, 1, 1, i ), 'linewidth', 4, 'color', 'r' );
    qy{ i } = quiver3( a, p_data_arr( 1, 1, i ), p_data_arr( 2, 1, i ), p_data_arr( 3, 1, i ), 2*tmp_scl * R_data_arr( 1, 2, 1, i ), 2*tmp_scl * R_data_arr( 2, 2, 1, i ), 2*tmp_scl * R_data_arr( 3, 2, 1, i ), 'linewidth', 4, 'color', 'g' );
    qz{ i } = quiver3( a, p_data_arr( 1, 1, i ), p_data_arr( 2, 1, i ), p_data_arr( 3, 1, i ), 2*tmp_scl * R_data_arr( 1, 3, 1, i ), 2*tmp_scl * R_data_arr( 2, 3, 1, i ), 2*tmp_scl * R_data_arr( 3, 3, 1, i ), 'linewidth', 4, 'color', 'b' );
end

% Draw the main marker
sc = scatter3( a, xc_arr( 2, 1 ), xc_arr( 3, 1 ), xc_arr( 4, 1 ), 100, 'filled', 'o', 'markerfacecolor', 'white', 'markeredgecolor', 'black', 'linewidth', 4 );

for i = 1 : Nstep : length( t_arr )
    
    for j = 1 : Ntraj
        set( s_arr{ j }, 'XData', p_data_arr( 1, i, j ), 'YData', p_data_arr( 2, i, j ), 'ZData', p_data_arr( 3, i, j ) );
        set( qx{ j }, 'XData', p_data_arr( 1, i, j ), 'YData', p_data_arr( 2, i, j ), 'ZData', p_data_arr( 3, i, j ), 'UData', 2*tmp_scl * R_data_arr( 1, 1, i, j ), 'VData', 2*tmp_scl * R_data_arr( 2, 1, i, j ), 'WData', 2*tmp_scl * R_data_arr( 3, 1, i, j ))
        set( qz{ j }, 'XData', p_data_arr( 1, i, j ), 'YData', p_data_arr( 2, i, j ), 'ZData', p_data_arr( 3, i, j ), 'UData', 2*tmp_scl * R_data_arr( 1, 2, i, j ), 'VData', 2*tmp_scl * R_data_arr( 2, 2, i, j ), 'WData', 2*tmp_scl * R_data_arr( 3, 2, i, j ))
        set( qy{ j }, 'XData', p_data_arr( 1, i, j ), 'YData', p_data_arr( 2, i, j ), 'ZData', p_data_arr( 3, i, j ), 'UData', 2*tmp_scl * R_data_arr( 1, 3, i, j ), 'VData', 2*tmp_scl * R_data_arr( 2, 3, i, j ), 'WData', 2*tmp_scl * R_data_arr( 3, 3, i, j ))
    end
    
    set( sc, 'XData', xc_arr( 2, i ), 'YData', xc_arr( 3, i ), 'ZData', xc_arr( 4, i ) )
   

    drawnow 
    
    set( t1, 'string', sprintf( 'Time %.3f s', t_arr( i ) ) );
    
    tmp_frame = getframe( f );
%     writeVideo( v,tmp_frame );
    i
end
% close( v );

%%  -- (1B) Generating the trajectories for Position 

