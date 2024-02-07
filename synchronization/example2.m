% [Title]     Study of Synchronization, Using Complex Number5
% [Author]    Moses Chong-ook Nah
% [Email]     mosesnah@mit.edu
% [Update]    At 2024.02.02

%% (--) Initialization
clear; close all; clc;
fig_config( 'fontSize', 20, 'markerSize', 10 )

%% (1A) Plot the Andronov Hopf Oscillator 
close all;

% Parameters of the Andronov-Hopf Oscillator 
r = 1.0; w = 2*pi;
f = figure( ); a = axes( 'parent', f );
hold on; axis equal
N_tmp = 30;

% First, run the simulation for both oscillators
dt = 1e-3;
T  = 12.0;
N  = round( T/dt );

for j = 1 : N_tmp
    % Using complex number
    z = 2*complex( rand, rand )-complex( rand, rand ); 
    % Time array 
    t_arr = (0:N)*dt;
    z_arr = zeros( 1, N+1 );
    z_arr( 1 ) = z;
    for i = 2 : N+1
        dz = ( (r^2-abs(z)^2 ) + 1i*w ) * z;
        z = z + dz*dt;
    
        z_arr( i ) = z;
    end
    plot( a, real( z_arr ), imag( z_arr ), 'linewidth', 3, 'color', 'k'  );
end

lw = 3.0;
set( a, 'xlim', [-lw, lw], 'ylim', [-lw, lw ])

%% (1B) Andronov Hopf Oscillators with Two Limit Cycles

% Now, consider a Oscillator with two limit cycles
% We use three Andronov-Hopf Oscillator 
phi1 = zeros( 3, 1 );
phi2 = [ 0, 2*pi/3, 4*pi/3]';

% Get the X Matrix
X = exp( 1i*[phi1, phi2]);
N = null( X', 'r'); 

% Define a Hurwitz Matrix
V = -1;

% The M matrix
M = N * V * N';

% Now, define three complex numbers for the three oscillators 
% Parameters of the Andronov-Hopf Oscillator 
r = 1.0; w = 2*pi;

% Using complex number
r0 = r*0.1;
a  = 0.6;
phi = phi1 * a + phi2 * (1-a);
z_tmp = r0*exp( 1i*phi );
z1 = z_tmp( 1 ); 
z2 = z_tmp( 2 );
z3 = z_tmp( 3 );

% First, run the simulation for both oscillators
dt = 1e-3;
T  = 12.0;
N  = round( T/dt );

% Time array 
t_arr = (0:N)*dt;
z_arr = zeros( 3, N+1 );
z_arr( :, 1 ) = [ z1; z2; z3 ];

for i = 2 : N+1
    dz1 = ( (r^2-abs(z1)^2 ) + 1i*w ) * z1 + M( 1, : )*[z1;z2;z3];
    dz2 = ( (r^2-abs(z2)^2 ) + 1i*w ) * z2 + M( 2, : )*[z1;z2;z3];
    dz3 = ( (r^2-abs(z3)^2 ) + 1i*w ) * z3 + M( 3, : )*[z1;z2;z3];

    z1 = z1 + dz1*dt;
    z2 = z2 + dz2*dt;
    z3 = z3 + dz3*dt;

    z_arr( 1, i ) = z1;
    z_arr( 2, i ) = z2;
    z_arr( 3, i ) = z3;
end

f = figure( ); a = axes( 'parent', f );
hold on
plot( a, t_arr, real( z_arr( 1, : ) ) )
plot( a, t_arr, real( z_arr( 2, : ) ) )
plot( a, t_arr, real( z_arr( 3, : ) ) )

