% [Title]     Plotting the Transformation System
% [Author]    Moses Chong-ook Nah
% [Email]     mosesnah@mit.edu
% [Update]    At 2023.11.12

%% (--) Initialization
clear; close all; clc;
fig_config( 'fontSize', 20, 'markerSize', 10 )

%% =======================================================
%% (1-) Transformation System
%%  -- (1A) Transformation System for Discrete movement

f = figure( ); a = axes( 'parent', f );
hold( a, 'on' );

for tmp = [0.1,0.25, 0.5, 0.75, 1.0, 2.5]

    % Canonical System
    tau     = tmp;
    cs      = CanonicalSystem( 'discrete', tau, 1 );

    % Gains of the Transformation System
    alpha_z = 1.0;
    beta_z  = alpha_z/4;
    g       = 2.0;
    y0      = 1.0;
    z0      = 0.0;
    T       = 10;
    
    % The number of time steps for the intergration
    N     = 1000;
    t_arr = linspace( 0, T, N+1 ); 
    
    trans_sys = TransformationSystem( alpha_z, beta_z, cs );
    [ y_arr, z_arr, ~ ] = trans_sys.rollout( y0, z0, g, zeros( 1, N ), 1, t_arr );
    
    hold on
    plot( a, t_arr, y_arr, 'color',  [0.4660 0.6740 0.1880] )
    set( a, 'xticklabel', { }, 'yticklabel', { } )

end
set( a, 'ylim', [ 0.0, 3.0 ] )

%%  -- (1B) Transformation System for Rhythmic movement
% TBD