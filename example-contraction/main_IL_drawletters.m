% [Title]     Saving the Weights for Imitation Learning, Import DAta
% [Author]    Moses Chong-ook Nah
% [Email]     mosesnah@mit.edu
% [Update]    At 2023.11.15

%% (--) Initialization
clear; close all; clc;
fig_config( 'fontSize', 20, 'markerSize', 10 )

%% =====================================================================
%% (1-) Learning the weights for the Imitation Learning
%%  -- (1A) Import Data

% Import the learned weights
alphabets = { 'R', 'A', 'L' };
Na = length( alphabets );

% dataset 
traj_data = cell( 1, Na );

for i = 1 : Na
    a = alphabets{ i };

    % Call Data
    load( [ './learned_parameters/', a, '.mat' ] );
    traj_data{ i } = data;
end

%%  -- (1B) Draw the Letters

t0i   = 0.0;
T     = 8;
dt    = 1e-4;
t_arr = 0:dt:T;
Nt    = length( t_arr );

p_data_arr = zeros( 2, Nt, Na );

scl_arr = [ 1.0, 1.5, 1.3];

% We need to first generate the DMP 
for i = 1 : Na
    
    % Call the data from the trajectory
    data = traj_data{ i };
    
    scl = scl_arr( i );

    % Getting the number of basis functions from the nonlinear forcing term
    [ ~, N ] = size( data.weight );

    % The Three elements of DMP
    cs        = CanonicalSystem( 'discrete', data.tau, data.alpha_s );
    fs        = NonlinearForcingTerm( cs, N );
    trans_sys = TransformationSystem( data.alpha_z, data.beta_z, cs );

    % Calculate the nonlinear forcing term 
    % This can be diminishing
    input_arr = fs.calc_forcing_term( t_arr( 1:end-1 ), data.weight, t0i, eye( 2 ) );

    % Rollout, note that initial condition is set to be zeros. 
    [ y_arr, z_arr, dy_arr ] = trans_sys.rollout( zeros( 2, 1 ), zeros( 2, 1 ),scl *data.goal, scl * input_arr, t0i, t_arr  );    
    
     p_data_arr( :, :, i ) =  y_arr;

end

%% -- (1C) Plotting the alphabets

xy_off = [ 0,0 ; 10, -1; 23, 12 ]';

f = figure( ); a = axes( 'parent', f );
hold( a, 'on' )
for i = 1 : Na
    tmp = p_data_arr( :, :, i ) + xy_off( :, i );
    plot( a, tmp( 1, : ), tmp( 2, : ), 'linewidth', 4, 'color' ,'k')
end
axis equal

%% =====================================================================
%% (2-) Using Contraction Theory
%%  -- (2A) Call the DAta

% Import the learned weights
alphabets = { 'R', 'A', 'L' };
Ntraj = length( alphabets );

% dataset 
traj_data = cell( 1, Ntraj );

for i = 1 : Ntraj
    a = alphabets{ i };

    % Call Data
    load( [ './learned_parameters/', a, '.mat' ] );
    traj_data{ i } = data;

    % Change the Tau to faster movements
    traj_data{ i }.tau = 3;
end

c_arr = [ 0.4940, 0.1840, 0.5560;
          0.6350, 0.0780, 0.1840;
      	  0.0000, 0.4470, 0.7410];	


%%  -- (2B) Generating the trajectories 

tinit =    1.0;           % The initial time of the simulation
T     =   10.0;           % The   whole time of the simulation 
toff  =   -0.1;           % Time-offset for the subsequent movement to start, must be negative
dt    =   1e-3;           % Time-step  for the simulation
t_arr = 0:dt:T;           % Time array for the simulation

% Time offset must be strictly negative. 
assert( toff < 0 );

% Position and velocity array of the trajectory
% This is 2D Data
 p_data_arr = zeros( 2, length( t_arr ), Ntraj );
dp_data_arr = zeros( 2, length( t_arr ), Ntraj );

% Scaling and Rotation of the trajectory
scl_arr = [ 1.0, 1.5, 1.3];

% The input force array for the movements
force_arr = zeros( 2, length( t_arr )-1, Ntraj );

% The goal and initated time and final time for each movement
goal_arr = zeros( 2, Ntraj );
 t0i_arr = zeros( 1, Ntraj );
 t0f_arr = zeros( 1, Ntraj );

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
    fs        = NonlinearForcingTerm( cs, N );
    trans_sys = TransformationSystem( data.alpha_z, data.beta_z, cs );

    % Calculate the nonlinear forcing term 
    % This can be diminishing
    input_arr_trimmed = fs.calc_forcing_term( t_arr( 1:end-1 ), data.weight, t0i, eye( 2 ), 'trimmed' );
    input_arr         = fs.calc_forcing_term( t_arr( 1:end-1 ), data.weight, t0i, eye( 2 ) );

    % Scaling
    scl = scl_arr( i );
        
    % Calculate the new goal location and save
    new_goal = scl*data.goal;
    goal_arr( :, i ) = new_goal;

    % Calculate the new force and save
    new_force = scl*input_arr; 
    new_force_trimmed = scl*input_arr_trimmed; 

    if i == Ntraj
        force_arr( :, :, i ) = new_force_trimmed;
    else
        force_arr( :, :, i ) = new_force;
    end

    % Rollout, note that initial condition is set to be zeros. 
    [ y_arr, z_arr, dy_arr ] = trans_sys.rollout( zeros( 2, 1 ), zeros( 2, 1 ), new_goal, new_force_trimmed, t0i, t_arr  );    
    
     p_data_arr( :, :, i ) =  y_arr;
    dp_data_arr( :, :, i ) = dy_arr;

    % Define the next movement's starting time.
    t0i = t0i + toff + data.tau;

end

% Add position offset
xy_off = [ 0,0 ;9, 0; 23, 12.5 ]';

p_data_arr( :, :, 2 ) = p_data_arr( :, :, 2 ) + xy_off( :, 2 );
p_data_arr( :, :, 3 ) = p_data_arr( :, :, 3 ) + xy_off( :, 3 );

goal_arr( :, 1 ) = p_data_arr( :, end, 1 );
goal_arr( :, 2 ) = p_data_arr( :, end, 2 );
goal_arr( :, 3 ) = p_data_arr( :, end, 3 );

% Plotting the positions to double cehck
f = figure( ); a = axes( 'parent', f );
axis equal; 
hold( a, 'on' )

for i = 1 : Ntraj
    plot( a, p_data_arr( 1, :, i ), p_data_arr( 2, :, i ), 'linewidth', 5, 'color', c_arr( i, : ) )
    scatter( a, p_data_arr( 1,   1, i ), p_data_arr( 2,   1, i ), 300, 'filled', 'o', 'markerfacecolor', 'white', 'markeredgecolor', c_arr( i, : ), 'linewidth', 4 )
    scatter( a, p_data_arr( 1, end, i ), p_data_arr( 2, end, i ), 300, 'filled', 'd', 'markerfacecolor', 'white', 'markeredgecolor', c_arr( i, : ), 'linewidth', 4 )    
end

%%  -- (2C) Use Contraction Theory

% Setting the Departure and Arrival Time 
t_depart_arr = t0f_arr( 1:2 ) - [ 0.1, 0.08 ];
t_arrive_arr = t0i_arr( 2:3 ) + [ 0.1, 0.30 ];

% Define the A matrix and az, bz values
% All should be identical, 
% Tau can be arbitrarily chosen based on the proof
% Choosing one arbitrarily
tau = traj_data{ 1 }.tau;
as  = traj_data{ 1 }.alpha_s;
az  = traj_data{ 1 }.alpha_z;
bz  = traj_data{ 1 }.beta_z;

% A Matrix
Amat = 1/tau * [            -as,       zeros( 1, 2 ),  zeros( 1, 2 ); 
                  zeros( 2, 1 ),       zeros( 2, 2 ),       eye( 2 );
                  zeros( 2, 1 ), -az * bz * eye( 2 ), -az * eye( 2 )];

Cg = 1/tau * [ zeros( 3, 1 ); az * bz * goal_arr( :, 1 )];

% The initial condition is identical to the first trajectory x0 value.
x0 = [ 1; p_data_arr( :, 1, 1 ); dp_data_arr( :, 1, 1 )*traj_data{1}.tau];

% The coupled x array 
xc_arr = zeros( 5, length( t_arr ) );
xc_arr( :, 1 ) = x0;

xcurr = x0;
% Iterating through the x_arr

% Integrating over time. 
for i = 1 : length( t_arr )-1
    
    t = t_arr( i );

    if t <= tinit
        dx = zeros( 5, 1 );
    else
        % Add the force for the first movement
        dx  = Amat * xcurr + 1/traj_data{ 1 }.tau * [ zeros( 3, 1 ); force_arr( :, i, 1 ) ] + Cg;

        % Iterating through the movement with activation
        for j = 2 : Ntraj
            fp = force_arr( :, i, j-1 );  % Force for previous movement
            fn = force_arr( :, i, j   );  % Force for the next movement
            
            Dp = traj_data{ j-1 }.tau;
            Dn = traj_data{ j   }.tau;

            td   = t_depart_arr( j-1 );
            ta   = t_arrive_arr( j-1 );
        
            gp = 1/Dp * [ zeros( 3, 1 ); az * bz * goal_arr( :, j-1 ) ];
            gn = 1/Dn * [ zeros( 3, 1 ); az * bz * goal_arr( :, j   ) ];
            
            gain = clip_func( t, td, ta );
            dx  = dx +  gain * ( - 1/Dp * [ zeros( 3, 1 ); fp ] + 1/Dn * [ zeros( 3, 1 ); fn ] ) ...
                     +  gain * ( -gp + gn ); 
        end

    end


    xc_arr( :, i+1 ) = xc_arr( :, i ) + dx * dt;     
    xcurr = xc_arr( :, i+1 );

end

%%  -- (2D) Image for RA-L

% Plotting the positions to double cehck
f = figure( ); a = axes( 'parent', f );
axis equal; 
hold( a, 'on' )
set( a, 'visible', 'off' )

for i = 1 : Ntraj
    plot( a, p_data_arr( 1, :, i ), p_data_arr( 2, :, i ), 'linewidth', 15, 'color', c_arr( i, : ) )
    scatter( a, p_data_arr( 1,   1, i ), p_data_arr( 2,   1, i ), 400, 'filled', 'o', 'markerfacecolor', 'white', 'markeredgecolor', c_arr( i, : ), 'linewidth', 8 )
    scatter( a, p_data_arr( 1, end, i ), p_data_arr( 2, end, i ), 400, 'filled', 'd', 'markerfacecolor', 'white', 'markeredgecolor', c_arr( i, : ), 'linewidth', 8 )    
end

% Overlap the plot with the contracte trajectory 
plot( a, xc_arr( 2, : ), xc_arr( 3, : ), 'linewidth', 6, 'color', 'black' )

fig_save( f, './figures/image1')

%%  -- (2E) Generating the Video    

% v = VideoWriter( 'video.mp4','MPEG-4' );
% v.FrameRate = 30;

% open( v );f
Nstep = round( 1/dt / 30 );

f = figure( ); a = axes( 'parent', f );
t1 = title( sprintf( 'Time %.3f s', 0 ) );
axis equal;
set( a, 'xlim', [ -5, 33 ], 'ylim', [-3, 15] )
hold on
s_arr = cell( 1, Ntraj );

for i = 1: Ntraj
    s_arr{ i } = scatter3( a, p_data_arr( 1, 1, i ), p_data_arr( 2, 1, i ), 100, 'filled', 'o', 'markerfacecolor', c_arr( i, : ), 'markeredgecolor', c_arr( i, : ), 'linewidth', 4 );

    plot( a, p_data_arr( 1, :, i ), p_data_arr( 2, :, i ), 'linewidth', 5, 'color', c_arr( i, : ) )
    scatter( a, p_data_arr( 1,   1, i ), p_data_arr( 2,   1, i ), 20, 'filled', 'o', 'markerfacecolor', 'white', 'markeredgecolor', c_arr( i, : ), 'linewidth', 4 )
    scatter( a, p_data_arr( 1, end, i ), p_data_arr( 2, end, i ), 20, 'filled', 'd', 'markerfacecolor', 'white', 'markeredgecolor', c_arr( i, : ), 'linewidth', 4 )  

end

% Draw the main marker
sc = scatter( a, xc_arr( 2, 1 ), xc_arr( 3, 1 ), 100, 'filled', 'o', 'markerfacecolor', 'white', 'markeredgecolor', 'black', 'linewidth', 4 );

for i = 1 : Nstep : length( t_arr )
    
    for j = 1 : Ntraj
        set( s_arr{ j }, 'XData', p_data_arr( 1, i, j ), 'YData', p_data_arr( 2, i, j ) );
    end
    
    set( sc, 'XData', xc_arr( 2, i ), 'YData', xc_arr( 3, i ) )

    drawnow 
    
    set( t1, 'string', sprintf( 'Time %.3f s', t_arr( i ) ) );
    
    tmp_frame = getframe( f );
%     writeVideo( v,tmp_frame );
    i
end
% close( v );

%% =====================================================================
%% (3-) Spectrum of Movements
%%  -- (3A) Call the Data

% Change from R and A to L

% Import the learned weights
alphabets = { 'R', 'A', 'L' };
Na = length( alphabets );

% dataset 
traj_data = cell( 1, Na );

for i = 1 : Na
    a = alphabets{ i };

    % Call Data
    load( [ './learned_parameters/', a, '.mat' ] );
    data.tau = 3;
    traj_data{ i } = data;
end


%%  -- (3B) Draw the Letters

t0i   = 0.0;
T     = 4;
dt    = 1e-4;
t_arr = 0:dt:T;
Nt    = length( t_arr );

p_data_arr = zeros( 2, Nt, Na );
force_arr  = zeros( 2, Nt-1, Na );
goal_arr   = zeros( 2, Na );
scl_arr    = [ 1.0, 1.5, 1.3];

% We need to first generate the DMP 
for i = 1 : Na
    
    % Call the data from the trajectory
    data = traj_data{ i };
    
    scl = scl_arr( i );

    % Getting the number of basis functions from the nonlinear forcing term
    [ ~, N ] = size( data.weight );

    % The Three elements of DMP
    cs        = CanonicalSystem( 'discrete', data.tau, data.alpha_s );
    fs        = NonlinearForcingTerm( cs, N );
    trans_sys = TransformationSystem( data.alpha_z, data.beta_z, cs );

    % Calculate the nonlinear forcing term 
    % This can be diminishing
    input_arr  = fs.calc_forcing_term( t_arr( 1:end-1 ), data.weight, t0i, eye( 2 ), 'trimmed' );
    input_arr2 = fs.calc_forcing_term( t_arr( 1:end-1 ), data.weight, t0i, eye( 2 ));

    % Rollout, note that initial condition is set to be zeros. 
    [ y_arr, z_arr, dy_arr ] = trans_sys.rollout( zeros( 2, 1 ), zeros( 2, 1 ),scl*data.goal, scl*input_arr, t0i, t_arr  );    
    
     p_data_arr(  :, :, i ) = y_arr;
     force_arr( :, :, i ) = scl*input_arr;
     goal_arr( :, i )     = scl*data.goal;

end

% For the third movement, we also match the y trajectory
p_data_arr(  :, :, 3 ) = p_data_arr(  :, :, 3 ) - p_data_arr(  :, end, 3 );

%%  -- (1C) Spectrum of Movement

% There are R, A and L
% Number of spectrum
Ns = 16;
w_arr = linspace( 0, 1, Ns );

% Generate the first half of the array (1 to 0.3)
firstHalf = linspace(1, 0.3, ceil(Ns/2));

% Generate the second half of the array (0.3 to 1)
secondHalf = linspace(0.3, 1, floor(Ns/2));

% Concatenate the two halves to form a symmetrical array
scl_arr = [firstHalf, secondHalf];

offset_arr = 9*scl_arr;
offset_arr = cumsum( offset_arr );

p_data_arr_spectrum = zeros( 2, Nt, Ns );

f = figure( );  a = axes( 'parent', f );
hold( a, 'on' ); set( a, 'visible', 'off' )
axis equal

c_arr = [ 0.4940, 0.1840, 0.5560;
          0.6350, 0.0780, 0.1840;
      	  0.0000, 0.4470, 0.7410];	

for i = 1 : Ns

    % The Three elements of DMP
    cs        = CanonicalSystem( 'discrete', data.tau, data.alpha_s );
    fs        = NonlinearForcingTerm( cs, N );
    trans_sys = TransformationSystem( data.alpha_z, data.beta_z, cs );

    w = w_arr( i );

    % Calculate the nonlinear forcing term 
    % This can be diminishing
    goal1 = goal_arr( :, 1 );
    goal2 = goal_arr( :, 2 );

    scl = scl_arr( i );

    new_goal = scl * ((1-w) * goal1 + w * goal2);

    c1 = c_arr( 1, : );
    c2 = c_arr( 2, : );

    force1 = force_arr( :, :, 1 );
    force2 = force_arr( :, :, 2 );
    new_force = scl*((1-w) * force1 + w * force2);

    % Rollout, note that initial condition is set to be zeros. 
    [ y_arr, z_arr, dy_arr ] = trans_sys.rollout( zeros( 2, 1 ), zeros( 2, 1 ), new_goal, new_force, 0, t_arr  );    
    
    new_color = c1*(1-w) + c2*w;
    if i == 1 || i == Ns
        lw = 6;
    else
        lw = 3;
    end

    plot( a, y_arr( 1,: )+offset_arr( i ), y_arr( 2,: ), 'color', new_color, 'linewidth', lw )

end

% From A to L
offset_arr  = offset_arr + offset_arr( end ) - offset_arr( 1 );
offset_arr(2:end) =offset_arr(2:end)+ 2; 
for i = 1 : Ns

    % The Three elements of DMP
    cs        = CanonicalSystem( 'discrete', data.tau, data.alpha_s );
    fs        = NonlinearForcingTerm( cs, N );
    trans_sys = TransformationSystem( data.alpha_z, data.beta_z, cs );

    w = w_arr( i );

    % Calculate the nonlinear forcing term 
    % This can be diminishing
    goal2 = goal_arr( :, 2 );
    goal3 = goal_arr( :, 3 );

    scl = scl_arr( i );

    new_goal = scl * ((1-w) * goal2 + w * goal3);

    c1 = c_arr( 2, : );
    c2 = c_arr( 3, : );

    force1 = force_arr( :, :, 2 );
    force2 = force_arr( :, :, 3 );
    new_force = scl*((1-w) * force1 + w * force2);

    % Rollout, note that initial condition is set to be zeros. 
    [ y_arr, z_arr, dy_arr ] = trans_sys.rollout( zeros( 2, 1 ), zeros( 2, 1 ), new_goal, new_force, 0, t_arr  );    
    
    new_color = (1-w)*c1+w*c2;
        
    if i == 1 
        tmp = 0;
    else
        tmp = -min( y_arr( 2,: ) );
    end

    if i == 1 || i == Ns
        lw = 6;
    else
        lw = 3;
    end
    
    plot( a, y_arr( 1,: )+offset_arr( i ), y_arr( 2,: )+tmp, 'color', new_color, 'linewidth', lw )
    
end

fig_save( f, './figures/image2')

