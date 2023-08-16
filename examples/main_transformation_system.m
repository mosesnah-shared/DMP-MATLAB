%% [Example Script] The basic Transformation System
% [Author] Moses Chong-ook Nah
% [Date]   2023.08.15

%% [Initialization] 
clear; close all; clc;
fig_config( 'fontSize', 20, 'markerSize', 10 )

%% [Transformation System]
% The basic transformation system which is a second-order linear system

alpha_z = 1.0;
beta_z  = 1.0;
tau     = 1.0;
g       = 1.0;
y0      = 0.0;
z0      = 0.0;

trans_sys = TransformationSystem( alpha_z, beta_z, tau, y0, z0 );

dt = 0.001;
Nt = 50000;

t_arr = dt * 0:Nt;
y_arr = zeros( 1, Nt + 1 );
z_arr = zeros( 1, Nt + 1 );

y_arr( 1 ) = y0;
z_arr( 1 ) = z0;

for i = 2 : Nt+1
    [ y, z, ~, ~ ] = trans_sys.step( g, 0, dt );
    y_arr( i ) = y;
    z_arr( i ) = z;
end

hold on
plot( t_arr, y_arr )
plot( t_arr, z_arr )