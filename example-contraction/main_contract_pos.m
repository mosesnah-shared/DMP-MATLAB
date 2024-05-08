% [Title]     Main Contraction Theory Position Example
% [Author]    Moses Chong-ook Nah
% [Email]     mosesnah@mit.edu
% [Update]    At 2023.11.16

%% (--) Initialization
clear; close all; clc;
fig_config( 'fontSize', 20, 'markerSize', 10 )

%% =======================================================
%% (1-) Data Collection and Trajectory Generation
%%  -- (1A) Calling the trajectory

% The string list for the types of trajectories
% We have four trajectories for this example
traj_names = [ "cosine", "radial", "drawM_pos", "lift_up_down_pos" ];
Ntraj = length( traj_names );

traj_data = cell( 1, Ntraj );

% Loading the data
for i = 1 : Ntraj
    tmp = load( [ './learned_parameters/discrete/', traj_names{ i }, '.mat' ] );
    traj_data{ i } = tmp.data;
end

% We need to check this
traj_data{3}.tau = 3;
traj_data{4}.tau = 3;

% Color array for plots
c_arr = [ 0.0000, 0.4470, 0.7410;
          0.8500, 0.3250, 0.0980;
          0.9290, 0.6940, 0.1250;
          0.4940, 0.1840, 0.5560;
          0.4660, 0.6740, 0.1880;
          0.3010, 0.7450, 0.9330];

clear tmp*

%%  -- (1B) Generating the trajectories 

tinit =    1.0;           % The initial time of the simulation
T     =   15.0;           % The   whole time of the simulation 
toff  =   -0.5;           % Time-offset for the subsequent movement to start, must be negative
dt    =   1e-3;           % Time-step  for the simulation
t_arr = 0:dt:T;           % Time array for the simulation

% Time offset must be strictly negative. 
assert( toff < 0 );

% Position and velocity array of the trajectory
 p_data_arr = zeros( 3, length( t_arr ), Ntraj );
dp_data_arr = zeros( 3, length( t_arr ), Ntraj );

% Scaling and Rotation of the trajectory
scl_arr = [ 0.7, 1.0, 2.5, 5.0 ];

R_arr = zeros( 3, 3, 4 );
R_arr( :, :, 1 ) =  rotx(  20 ) * roty( -30 ) * rotz(  10 );
R_arr( :, :, 2 ) =  rotx(  10 ) * roty(  40 ) * rotz( -45 );
R_arr( :, :, 3 ) =  rotx(  45 ) * roty( -30 ) * rotz(  40 );
R_arr( :, :, 4 ) =  rotx( -45 ) * roty(  45 ) * rotz( -45 );

% The input force array for the movements
force_arr = zeros( 3, length( t_arr )-1, 4 );

% The goal and initated time and final time for each movement
goal_arr = zeros( 3, 4 );
 t0i_arr = zeros( 1, 4 );
 t0f_arr = zeros( 1, 4 );

% The first movement starts at tinit, hence setting up
 t0i = tinit;          

% We need to first generate the DMP 
for i = 1 : Ntraj
    
    % Call the data from the trajectory
    data = traj_data{ i };

    % Saving the initial time in the array
    t0i_arr( i ) = t0i;
    t0f_arr( i ) = t0i + data.tau;
    
    % Getting the number of basis functions from the nonlinear forcing term
    [ ~, N ] = size( data.weight );

    % The Three elements of DMP
    cs        = CanonicalSystem( 'discrete', data.tau, data.alpha_s );
    fs        =     gTerm( cs, N );
    trans_sys = TransformationSystem( data.alpha_z, data.beta_z, cs );

    % Calculate the nonlinear forcing term 
    % This can be diminishing
    input_arr = fs.calc_forcing_term( t_arr( 1:end-1 ), data.weight, t0i, eye( 3 ) );

    % Rotation and Scaling
     R  = R_arr( :, :, i );
    scl = scl_arr( i );
        
    % Calculate the new goal location and save
    new_goal = scl*R*data.goal;
    goal_arr( :, i ) = new_goal;

    % Calculate the new force and save
    new_force = scl*R*input_arr; 
    force_arr( :, :, i ) = new_force;

    % Rollout, note that initial condition is set to be zeros. 
    [ y_arr, z_arr, dy_arr ] = trans_sys.rollout( zeros( 3, 1 ), zeros( 3, 1 ), new_goal, new_force, t0i, t_arr  );    
    
     p_data_arr( :, :, i ) =  y_arr;
    dp_data_arr( :, :, i ) = dy_arr;

    % Define the next movement's starting time.
    t0i = t0i + toff + data.tau;

end

% Add position offset 
poff = zeros( 3, 3 );
poff( :, 1 ) = [  0.1; 0.2;  0.1 ];
poff( :, 2 ) = [ -0.1; 0.2; -0.3 ];
poff( :, 3 ) = [  0.1; 0.1;  0.0 ];

p_data_arr( :, :, 2 ) = p_data_arr( :, :, 2 ) + goal_arr( :, 1 ) + poff( :, 1 );
p_data_arr( :, :, 3 ) = p_data_arr( :, :, 3 ) + goal_arr( :, 2 ) + poff( :, 2 );
p_data_arr( :, :, 4 ) = p_data_arr( :, :, 4 ) + goal_arr( :, 3 ) + poff( :, 3 );

goal_arr( :, 1 ) = p_data_arr( :, end, 1 );
goal_arr( :, 2 ) = p_data_arr( :, end, 2 );
goal_arr( :, 3 ) = p_data_arr( :, end, 3 );
goal_arr( :, 4 ) = p_data_arr( :, end, 4 );


% Plotting the positions to double cehck
f = figure( ); a = axes( 'parent', f );
axis equal; view( 3 );
hold( a, 'on' )

for i = 1 : Ntraj
    plot3( a, p_data_arr( 1, :, i ), p_data_arr( 2, :, i ), p_data_arr( 3, :, i ), 'linewidth', 5, 'color', c_arr( i, : ) )
    scatter3( a, p_data_arr( 1,   1, i ), p_data_arr( 2,   1, i ), p_data_arr( 3,   1, i ), 'filled', 'o', 'markerfacecolor', 'white', 'markeredgecolor', c_arr( i, : ), 'linewidth', 4 )
    scatter3( a, p_data_arr( 1, end, i ), p_data_arr( 2, end, i ), p_data_arr( 3, end, i ), 'filled', 'd', 'markerfacecolor', 'white', 'markeredgecolor', c_arr( i, : ), 'linewidth', 4 )    
end

%%  -- (1C) Use Contraction Theory

% Setting the Departure and Arrival Time 
t_depart_arr = t0f_arr( 1:3 ) + toff + 0.1;
t_arrive_arr = t0i_arr( 2:4 ) + 0.5;

% Define the A matrix and az, bz values
% All should be identical, 
% Tau can be arbitrarily chosen based on the proof
% Choosing one arbitrarily
tau = traj_data{ 1 }.tau;
as  = traj_data{ 1 }.alpha_s;
az  = traj_data{ 1 }.alpha_z;
bz  = traj_data{ 1 }.beta_z;

% A Matrix
Amat = 1/tau * [            -as,       zeros( 1, 3 ),  zeros( 1, 3 ); 
                  zeros( 3, 1 ),       zeros( 3, 3 ),       eye( 3 );
                  zeros( 3, 1 ), -az * bz * eye( 3 ), -az * eye( 3 )];

Cg = 1/tau * [ zeros( 4, 1 ); az * bz * goal_arr( :, 1 )];

% The initial condition is identical to the first trajectory x0 value.
x0 = [ 1; p_data_arr( :, 1, 1 ); dp_data_arr( :, 1, 1 )*traj_data{1}.tau];

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
        dx  = Amat * xcurr + 1/traj_data{ 1 }.tau * [ zeros( 4, 1 ); force_arr( :, i, 1 ) ] + Cg;

        % Iterating through the movement with activation
        for j = 2 : 4
            fp = force_arr( :, i, j-1 );  % Force for previous movement
            fn = force_arr( :, i, j   );  % Force for the next movement
            
            Dp = traj_data{ j-1 }.tau;
            Dn = traj_data{ j   }.tau;

            td   = t_depart_arr( j-1 );
            ta   = t_arrive_arr( j-1 );
        
            gp = 1/Dp * [ zeros( 4, 1 ); az * bz * goal_arr( :, j-1 ) ];
            gn = 1/Dn * [ zeros( 4, 1 ); az * bz * goal_arr( :, j   ) ];
            
            gain = clip_func( t, td, ta );
            dx  = dx +  gain * ( - 1/Dp * [ zeros( 4, 1 ); fp ] + 1/Dn * [ zeros( 4, 1 ); fn ] ) ...
                     +  gain * ( -gp + gn ); 
        end

    end


    xc_arr( :, i+1 ) = xc_arr( :, i ) + dx * dt;     
    xcurr = xc_arr( :, i+1 );

end


%%  -- (1E) Generating the Video    

% v = VideoWriter( 'video.mp4','MPEG-4' );
% v.FrameRate = 30;

% open( v );
Nstep = round( 1/dt / 30 );

f = figure( ); a = axes( 'parent', f );
t1 = title( sprintf( 'Time %.3f s', 0 ) );

axis equal; view( 3 );
hold on
s_arr = cell( 1, Ntraj );

for i = 1: Ntraj
    s_arr{ i } = scatter3( a, p_data_arr( 1, 1, i ), p_data_arr( 2, 1, i ), p_data_arr( 3, 1, i ), 100, 'filled', 'o', 'markerfacecolor', c_arr( i, : ), 'markeredgecolor', c_arr( i, : ), 'linewidth', 4 );

    plot3( a, p_data_arr( 1, :, i ), p_data_arr( 2, :, i ), p_data_arr( 3, :, i ), 'linewidth', 5, 'color', c_arr( i, : ) )
    scatter3( a, p_data_arr( 1,   1, i ), p_data_arr( 2,   1, i ), p_data_arr( 3,   1, i ), 20, 'filled', 'o', 'markerfacecolor', 'white', 'markeredgecolor', c_arr( i, : ), 'linewidth', 4 )
    scatter3( a, p_data_arr( 1, end, i ), p_data_arr( 2, end, i ), p_data_arr( 3, end, i ), 20, 'filled', 'd', 'markerfacecolor', 'white', 'markeredgecolor', c_arr( i, : ), 'linewidth', 4 )  

end

% Draw the main marker
sc = scatter3( a, xc_arr( 2, 1 ), xc_arr( 3, 1 ), xc_arr( 4, 1 ), 100, 'filled', 'o', 'markerfacecolor', 'white', 'markeredgecolor', 'black', 'linewidth', 4 );

for i = 1 : Nstep : length( t_arr )
    
    for j = 1 : Ntraj
        set( s_arr{ j }, 'XData', p_data_arr( 1, i, j ), 'YData', p_data_arr( 2, i, j ), 'ZData', p_data_arr( 3, i, j ) );
    end
    
    set( sc, 'XData', xc_arr( 2, i ), 'YData', xc_arr( 3, i ), 'ZData', xc_arr( 4, i ) )

    drawnow 
    
    set( t1, 'string', sprintf( 'Time %.3f s', t_arr( i ) ) );
    
    tmp_frame = getframe( f );
%     writeVideo( v,tmp_frame );
    i
end
% close( v );
