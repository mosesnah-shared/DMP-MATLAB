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
traj_names = [ "cosine", "radial", "drawM_pos", "lift_up_down_pos" ];
Ntraj = length( traj_names );

traj_data = cell( 1, Ntraj );

for i = 1 : Ntraj
    tmp = load( [ './learned_parameters/', traj_names{ i }, '.mat' ] );
    traj_data{ i } = tmp.data;
end

c_arr = [ 0.0000, 0.4470, 0.7410;
          0.8500, 0.3250, 0.0980;
          0.9290, 0.6940, 0.1250;
          0.4940, 0.1840, 0.5560;
          0.4660, 0.6740, 0.1880;
          0.3010, 0.7450, 0.9330];

%%  -- (1B) Generating the trajectories for Check

tinit =  1.0;
t0i   = tinit;
toff  =  -0.5;
T     = 25.0;
dt    =  1e-3;
t_arr = 0:dt:T;

 p_data_arr = zeros( 3, length( t_arr ), Ntraj );
dp_data_arr = zeros( 3, length( t_arr ), Ntraj );


scl_arr = [ 0.7, 1.0, 2.5, 5.0 ];

R_arr = zeros( 3, 3, 4 );
R_arr( :, :, 1 ) =  rotx(  20 ) * roty( -30 ) * rotz(  10 );
R_arr( :, :, 2 ) =  rotx(  10 ) * roty(  40 ) * rotz( -45 );
R_arr( :, :, 3 ) =  rotx(  45 ) * roty( -30 ) * rotz(  40 );
R_arr( :, :, 4 ) =  rotx( -45 ) * roty(  45 ) * rotz( -45 );

goal_arr = zeros( 3, 4 );
t0i_arr  = zeros( 1, 4 );

force_arr = zeros( 3, length( t_arr )-1, 4 );


% We need to first generate the DMP 
for i = 1 : Ntraj
    
    % Saving the initial time
    t0i_arr( i ) = t0i;

    data = traj_data{ i };
    
    % Getting the number of basis functions from the nonlinear forcing term
    [ ~, N ] = size( data.weight );

    % The Three elements of DMP
    cs        = CanonicalSystem( 'discrete', data.tau, data.alpha_s );
    fs        = NonlinearForcingTerm( cs, N );
    trans_sys = TransformationSystem( data.alpha_z, data.beta_z, cs );

    % The nonlinear forcing term array
    input_arr = fs.calc_forcing_term( t_arr( 1:end-1 ), data.weight, t0i, eye( 3 ) );

     R  = R_arr( :, :, i );
    scl = scl_arr( i );
        
    new_goal = scl*R*data.goal;

    [ y_arr, z_arr, dy_arr ] = trans_sys.rollout( zeros( 3, 1 ), zeros( 3, 1 ), new_goal, scl*R*input_arr, t0i, t_arr  );    
    
     p_data_arr( :, :, i ) =  y_arr;
    dp_data_arr( :, :, i ) = dy_arr;

    t0i = t0i + toff + data.tau;

    goal_arr( :, i ) = new_goal;
    force_arr( :, :, i ) = scl*R*input_arr; 
end

% Adding position offset 
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


% Plotting the arrays 
f = figure( ); a = axes( 'parent', f );
axis equal; view( 3 );
hold( a, 'on' )

for i = 1 : Ntraj
    plot3( a, p_data_arr( 1, :, i ), p_data_arr( 2, :, i ), p_data_arr( 3, :, i ), 'linewidth', 5, 'color', c_arr( i, : ) )
    scatter3( a, p_data_arr( 1,   1, i ), p_data_arr( 2,   1, i ), p_data_arr( 3,   1, i ), 'filled', 'o', 'markerfacecolor', 'white', 'markeredgecolor', c_arr( i, : ), 'linewidth', 4 )
    scatter3( a, p_data_arr( 1, end, i ), p_data_arr( 2, end, i ), p_data_arr( 3, end, i ), 'filled', 'd', 'markerfacecolor', 'white', 'markeredgecolor', c_arr( i, : ), 'linewidth', 4 )    
end

%%  -- (1C) Generating the state-space representation 

% Defining the four trajectories in state-space representation 
x_arr = zeros( 7, length( t_arr ), 4 );

for i = 1 : Ntraj
    
    data = traj_data{ i };
    
    % The Three elements of DMP
    cs = CanonicalSystem( 'discrete', data.tau, data.alpha_s );

    % First, we fill in the p_data_arr 
    x_arr( 2:4, :, i ) =  p_data_arr( :, :, i );
    x_arr( 5:7, :, i ) = dp_data_arr( :, :, i )*data.tau;

    % For the canonical system, we account for the initial time
    % Simply use t_arr, but with offset time 
    tmp_tarr = t_arr - t0i_arr( i ); 
    tmp_tarr( tmp_tarr <= 0 ) = 0;

    x_arr( 1, :, i ) = cs.calc( tmp_tarr );
end

%%  -- (1D) Use Contraction Theory

% Since the trajectory data is defined, we integrate the equation
% Define the departure time.
t_depart = t0i_arr( 1:3 ) + 0.2;

% Gain of the departure 
k_gain = 10 * ones( 1,3 );

% Define the A matrix and az, bz values
% All should be identical, 
% Tau can be arbitrarily chosen based on the proof
% Choosing one arbitrarily
tau = traj_data{ 1 }.tau;
as  = traj_data{ 1 }.alpha_s;
az  = traj_data{ 1 }.alpha_z;
bz  = traj_data{ 1 }.beta_z;

Amat = 1/tau * [            -as,       zeros( 1, 3 ),  zeros( 1, 3 ); 
                  zeros( 3, 1 ),       zeros( 3, 3 ),       eye( 3 );
                  zeros( 3, 1 ), -az * bz * eye( 3 ), -az * eye( 3 )];

Cg = 1/tau * [ zeros( 4, 1 ); az * bz * goal_arr( :, 1 )];

% The initial condition is identical to the first trajectory x0 value.
x0 = x_arr( :, 1, 1 );

% The coupled x array 
xc_arr = zeros( 7, length( t_arr ) );
xc_arr( :, 1 ) = x0;

xcurr = x0;
% Iterating through the x_arr
act1 = 0;
act2 = 0;
act3 = 0;

for i = 1 : length( t_arr )-1

    % Wait until the first 
    t = t_arr( i );
    if t <= tinit 
        dx = zeros( 7, 1 );
    else
        dx = Amat * xc_arr( :, i ) + 1/tau*[ zeros( 4, 1 ); force_arr( :, i, 1 ) ] + Cg;

        if t >= 3.5 && t <= 5.5
            act1 = act1 + 0.005;
            if act1 >= 1.0
                act1 = 1.0;
            end
            dx = dx + act1 * ( -  1/tau*[ zeros( 4, 1 ); force_arr( :, i, 1 ) ] +  1/tau*[ zeros( 4, 1 ); force_arr( :, i, 2 ) ] - 1/tau * [ zeros( 4, 1 ); az * bz * goal_arr( :, 1 )]  + 1/tau * [ zeros( 4, 1 ); az * bz * goal_arr( :, 2 )]  ) ;


        elseif t >= 5.5 && t <= 9.0
            act2 = act2 + 0.005;
            if act2 >= 1.0
                act2 = 1.0;
            end
            dx = dx + act2 * (  1/tau*[ zeros( 4, 1 ); force_arr( :, i, 3 ) ] - 1/tau * [ zeros( 4, 1 ); az * bz * goal_arr( :, 1 )] + 1/tau * [ zeros( 4, 1 ); az * bz * goal_arr( :, 3 )] + 10 * ( -xc_arr( :, i ) + x_arr( :, i, 3 ) ) ) ;
     
        elseif t >= 9.0
            act3 = act3 + 0.005;
            if act3 >= 1.0
                act3 = 1.0;
            end
            dx = dx + act3 * (  1/tau*[ zeros( 4, 1 ); force_arr( :, i, 4 ) ] - 1/tau * [ zeros( 4, 1 ); az * bz * goal_arr( :, 1 )] + 1/tau * [ zeros( 4, 1 ); az * bz * goal_arr( :, 4 )] + 10 * ( -xc_arr( :, i ) + x_arr( :, i, 3 ) ) ) ;
          
        end


        
    end

    xc_arr( :, i+1 ) = xc_arr( :, i ) + dx * dt;     

end


%%  -- (1E) Generating the Video    

% v = VideoWriter( 'video.mp4','MPEG-4' );
% v.FrameRate = 30;

% open( v );f
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
