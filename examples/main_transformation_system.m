%% ==================================================================
%% [Title] Example Script - Plotting the Transformation System, without Input
% Author: Moses Chong-ook Nah
%  Email: mosesnah@mit.edu
%   Date: 2023.08.15
%% ==================================================================

%% [0A] Initialization
clear; close all; clc;
fig_config( 'fontSize', 20, 'markerSize', 10 )

%% [1-] Transformation System

for tmp = [0.1,0.25, 0.5, 0.75, 1.0, 2.5]

    alpha_z = 1.0;
    beta_z  = alpha_z/4;
    tau     = tmp;
    cs      = CanonicalSystem( 'discrete', tau, 1 );
    g       = 2.0;
    y0      = 1.0;
    z0      = 0.0;
    T       = 10;
    
    N = 1000;
    t_arr = linspace( 0, T, N+1 ); 
    
    trans_sys = TransformationSystem( alpha_z, beta_z, cs );
    [ y_arr, z_arr, ~ ] = trans_sys.rollout( y0, z0, g, zeros( 1, N ), 1, t_arr );
    
    hold on
    plot( t_arr, y_arr, 'color',  [0.4660 0.6740 0.1880] )
    set( gca, 'xticklabel', { }, 'yticklabel', { } )

end
set( gca, 'ylim', [ 0.0, 3.0 ] )