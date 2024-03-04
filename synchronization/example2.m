% [Title]     Study of Synchronization, Using Complex Number5
% [Author]    Moses Chong-ook Nah
% [Email]     mosesnah@mit.edu
% [Update]    At 2024.02.02

%% (--) Initialization
clear; close all; clc;
fig_config( 'fontSize', 20, 'markerSize', 10 )

c_arr = cell( 1, 7 );
Nc = 7;
c_arr{ 1 } = "#0072BD";
c_arr{ 2 } = "#D95319";
c_arr{ 3 } = "#EDB120";
c_arr{ 4 } = "#7E2F8E";
c_arr{ 5 } = "#77AC30";
c_arr{ 6 } = "#4DBEEE";
c_arr{ 7 } = "#4DBEEE";

%% (1A) Plot + Video of the Andronov Hopf Oscillator 
close all;

% Parameters of the Andronov-Hopf Oscillator 
r = 1.0; w = 2*pi;

% First, run the simulation for both oscillators
dt = 1e-3;
T  = 4.0;
N  = round( T/dt );

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

% Generating the video from these values 
f = figure( ); 
a1 = subplot( 2, 1, 1 ); set( a1, 'parent', f );
hold on

plot( a1, real( z_arr ), imag( z_arr ), 'linewidth', 5, 'color', c_arr{1} )
p1 = scatter( a1, real( z_arr( 1 ) ), imag( z_arr( 1 ) ), 400, 'o', 'filled', 'linewidth', 5, 'markerfacecolor', 'w', 'markeredgecolor', c_arr{1} );
axis equal
lw = 3.0;
set( a1, 'xlim', [-lw, lw], 'ylim', [-lw, lw] )

a2 = subplot( 2, 1, 2 ); set( a2, 'parent', f );
hold on

plot( a2, t_arr, angle( z_arr ), 'linewidth', 5, 'color', c_arr{1} )
p2 = scatter( a2, t_arr( 1 ), angle( z_arr( 1 ) ), 400, 'o', 'filled', 'linewidth', 5, 'markerfacecolor', 'w', 'markeredgecolor', c_arr{1} );

% Set up the video writer
outputVideo = VideoWriter( 'videos/synchronization_imag/simple_osc.mp4', 'MPEG-4' );
outputVideo.FrameRate = 60;
open( outputVideo );

% Time per frame
numFrames    = round( T*outputVideo.FrameRate);
timePerFrame = T / numFrames;

for frameIdx = 1:numFrames

    % Current time for this frame
    currentTime = (frameIdx - 1) * timePerFrame;
    
    % Find the closest time in your timeArray (or interpolate if needed)
    [~,  idx ] = min( abs( t_arr - currentTime ) );

    set( p1, 'xdata',  real( z_arr( idx ) ), 'ydata', imag( z_arr( idx ) ) )
    set( p2, 'xdata',  t_arr( idx ), 'ydata', angle( z_arr( idx ) ) )

    % Capture frame
    frame = getframe( f );
    
    % Write frame to video
    writeVideo( outputVideo, frame );
end

% Close the video file
close( outputVideo );

% Optionally, close the figure
close( f );



%% (1B) Andronov Hopf Oscillators with Two Limit Cycles
close all;
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
M1 = zeros( 3, 3 );
M2 = N * V * N';
    
% Now, define three complex numbers for the three oscillators 
% Parameters of the Andronov-Hopf Oscillator 
r = 1.0; w = 2*pi;

% Using complex number
r0 = rand( 3,1 );
% a  = 0.2;
phi = 8*rand( 3,1 )-8;
z_tmp = r0.*exp( 1i*phi );
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
M = M1;

for i = 2 : N+1
    dz1 = ( (r^2-abs(z1)^2 ) + 1i*w ) * z1 + M( 1, : )*[z1;z2;z3];
    dz2 = ( (r^2-abs(z2)^2 ) + 1i*w ) * z2 + M( 2, : )*[z1;z2;z3];
    dz3 = ( (r^2-abs(z3)^2 ) + 1i*w ) * z3 + M( 3, : )*[z1;z2;z3];

    if i == round( N/3 )
        M=M2;
    end

    z1 = z1 + dz1*dt;
    z2 = z2 + dz2*dt;
    z3 = z3 + dz3*dt;

    z_arr( 1, i ) = z1;
    z_arr( 2, i ) = z2;
    z_arr( 3, i ) = z3;
end

% Do the subplot
f = figure( ); 
a1 = subplot( 2, 3, 1 ); set( a1, 'parent', f ); hold on; axis equal;
a2 = subplot( 2, 3, 2 ); set( a2, 'parent', f ); hold on; axis equal;
a3 = subplot( 2, 3, 3 ); set( a3, 'parent', f ); hold on; axis equal;
lw = 2;
set( a1, 'xlim', [-lw, lw], 'ylim', [-lw, lw] )
set( a2, 'xlim', [-lw, lw], 'ylim', [-lw, lw] )
set( a3, 'xlim', [-lw, lw], 'ylim', [-lw, lw] )

plot( a1, real( z_arr( 1, : ) ), imag( z_arr( 1, : ) ), 'linewidth', 5, 'color', c_arr{ 1 } );
plot( a2, real( z_arr( 2, : ) ), imag( z_arr( 2, : ) ), 'linewidth', 5, 'color', c_arr{ 2 } );
plot( a3, real( z_arr( 3, : ) ), imag( z_arr( 3, : ) ), 'linewidth', 5, 'color', c_arr{ 3 } );

p1 = scatter( a1, real( z_arr( 1, 1 ) ), imag( z_arr( 1, 1 ) ), 300, 'o', 'filled', 'markerfacecolor', 'w', 'markeredgecolor', c_arr{ 1 }, 'linewidth', 5 );
p2 = scatter( a2, real( z_arr( 2, 1 ) ), imag( z_arr( 2, 1 ) ), 300, 'o', 'filled', 'markerfacecolor', 'w', 'markeredgecolor', c_arr{ 2 }, 'linewidth', 5 );
p3 = scatter( a3, real( z_arr( 3, 1 ) ), imag( z_arr( 3, 1 ) ), 300, 'o', 'filled', 'markerfacecolor', 'w', 'markeredgecolor', c_arr{ 3 }, 'linewidth', 5 );


a4 = subplot( 2, 3, 4:6 ); set( a4, 'parent', f ); hold on; 
set( a4, 'xlim', [0, max(t_arr)] )
plot( a4, t_arr, angle( z_arr( 1, : ) ), 'linewidth', 5, 'color', c_arr{ 1 } );
plot( a4, t_arr, angle( z_arr( 2, : ) ), 'linewidth', 5, 'color', c_arr{ 2 } );
plot( a4, t_arr, angle( z_arr( 3, : ) ), 'linewidth', 5, 'color', c_arr{ 3 } );
p4 = scatter( a4, t_arr( 1 ), angle( z_arr( 1, 1 ) ), 300, 'o', 'filled', 'markerfacecolor', 'w', 'markeredgecolor', c_arr{ 1 }, 'linewidth', 5 );
p5 = scatter( a4, t_arr( 1 ), angle( z_arr( 2, 1 ) ), 300, 'o', 'filled', 'markerfacecolor', 'w', 'markeredgecolor', c_arr{ 2 }, 'linewidth', 5 );
p6 = scatter( a4, t_arr( 1 ), angle( z_arr( 3, 1 ) ), 300, 'o', 'filled', 'markerfacecolor', 'w', 'markeredgecolor', c_arr{ 3 }, 'linewidth', 5 );


% Set up the video writer
outputVideo = VideoWriter( 'videos/synchronization_imag/two_limit_cycles.mp4', 'MPEG-4' );
outputVideo.FrameRate = 60;
open( outputVideo );

% Time per frame
numFrames    = round( T*outputVideo.FrameRate);
timePerFrame = T / numFrames;

for frameIdx = 1:numFrames

    % Current time for this frame
    currentTime = (frameIdx - 1) * timePerFrame;
    
    % Find the closest time in your timeArray (or interpolate if needed)
    [~,  idx ] = min( abs( t_arr - currentTime ) );

    set( p1, 'xdata',  real( z_arr( 1, idx ) ), 'ydata', imag( z_arr( 1, idx ) ) )
    set( p2, 'xdata',  real( z_arr( 2, idx ) ), 'ydata', imag( z_arr( 2, idx ) ) )
    set( p3, 'xdata',  real( z_arr( 3, idx ) ), 'ydata', imag( z_arr( 3, idx ) ) )

    set( p4, 'xdata',  t_arr( idx ), 'ydata', angle( z_arr( 1, idx ) ) )
    set( p5, 'xdata',  t_arr( idx ), 'ydata', angle( z_arr( 2, idx ) ) )
    set( p6, 'xdata',  t_arr( idx ), 'ydata', angle( z_arr( 3, idx ) ) )


    % Capture frame
    frame = getframe( f );
    
    % Write frame to video
    writeVideo( outputVideo, frame );
end

% Close the video file
close( outputVideo );

% Optionally, close the figure
close( f );

%% (1C) Andronov Hopf Oscillators with One Limit Cycle
close all;
phi = [ 0, 2*pi/3, 4*pi/3]';

% Get the X Matrix
X = exp( 1i*phi);
N = null( X', 'r'); 

% Define a Hurwitz Matrix
V = -1;

% The M matrix
M1 = zeros( 3, 3 );
M2 = N * V * N';
    
% Now, define three complex numbers for the three oscillators 
% Parameters of the Andronov-Hopf Oscillator 
r = 1.0; w = 2*pi;

% Using complex number
r0 = rand( 3,1 );
phi = 8*rand( 3,1 )-8;
z_tmp = r0.*exp( 1i*phi );
z1 = z_tmp( 1 ); 
z2 = z_tmp( 2 );
z3 = z_tmp( 3 );

% First, run the simulation for both oscillators
dt = 1e-3;
T  = 8.0;
N  = round( T/dt );

% Time array 
t_arr = (0:N)*dt;
z_arr = zeros( 3, N+1 );
z_arr( :, 1 ) = [ z1; z2; z3 ];
M = M1;

for i = 2 : N+1
    dz1 = ( (r^2-abs(z1)^2 ) + 1i*w ) * z1 + M( 1, : )*[z1;z2;z3];
    dz2 = ( (r^2-abs(z2)^2 ) + 1i*w ) * z2 + M( 2, : )*[z1;z2;z3];
    dz3 = ( (r^2-abs(z3)^2 ) + 1i*w ) * z3 + M( 3, : )*[z1;z2;z3];

    if i == round( N/3 )
        M=M2;
    end

    z1 = z1 + dz1*dt;
    z2 = z2 + dz2*dt;
    z3 = z3 + dz3*dt;

    z_arr( 1, i ) = z1;
    z_arr( 2, i ) = z2;
    z_arr( 3, i ) = z3;
end

% Do the subplot
f = figure( ); 
a1 = subplot( 2, 3, 1 ); set( a1, 'parent', f ); hold on; axis equal;
a2 = subplot( 2, 3, 2 ); set( a2, 'parent', f ); hold on; axis equal;
a3 = subplot( 2, 3, 3 ); set( a3, 'parent', f ); hold on; axis equal;
lw = 2;
set( a1, 'xlim', [-lw, lw], 'ylim', [-lw, lw] )
set( a2, 'xlim', [-lw, lw], 'ylim', [-lw, lw] )
set( a3, 'xlim', [-lw, lw], 'ylim', [-lw, lw] )

plot( a1, real( z_arr( 1, : ) ), imag( z_arr( 1, : ) ), 'linewidth', 5, 'color', c_arr{ 1 } );
plot( a2, real( z_arr( 2, : ) ), imag( z_arr( 2, : ) ), 'linewidth', 5, 'color', c_arr{ 2 } );
plot( a3, real( z_arr( 3, : ) ), imag( z_arr( 3, : ) ), 'linewidth', 5, 'color', c_arr{ 3 } );

p1 = scatter( a1, real( z_arr( 1, 1 ) ), imag( z_arr( 1, 1 ) ), 300, 'o', 'filled', 'markerfacecolor', 'w', 'markeredgecolor', c_arr{ 1 }, 'linewidth', 5 );
p2 = scatter( a2, real( z_arr( 2, 1 ) ), imag( z_arr( 2, 1 ) ), 300, 'o', 'filled', 'markerfacecolor', 'w', 'markeredgecolor', c_arr{ 2 }, 'linewidth', 5 );
p3 = scatter( a3, real( z_arr( 3, 1 ) ), imag( z_arr( 3, 1 ) ), 300, 'o', 'filled', 'markerfacecolor', 'w', 'markeredgecolor', c_arr{ 3 }, 'linewidth', 5 );


a4 = subplot( 2, 3, 4:6 ); set( a4, 'parent', f ); hold on; 
set( a4, 'xlim', [0, max(t_arr)] )
plot( a4, t_arr, angle( z_arr( 1, : ) ), 'linewidth', 5, 'color', c_arr{ 1 } );
plot( a4, t_arr, angle( z_arr( 2, : ) ), 'linewidth', 5, 'color', c_arr{ 2 } );
plot( a4, t_arr, angle( z_arr( 3, : ) ), 'linewidth', 5, 'color', c_arr{ 3 } );
p4 = scatter( a4, t_arr( 1 ), angle( z_arr( 1, 1 ) ), 300, 'o', 'filled', 'markerfacecolor', 'w', 'markeredgecolor', c_arr{ 1 }, 'linewidth', 5 );
p5 = scatter( a4, t_arr( 1 ), angle( z_arr( 2, 1 ) ), 300, 'o', 'filled', 'markerfacecolor', 'w', 'markeredgecolor', c_arr{ 2 }, 'linewidth', 5 );
p6 = scatter( a4, t_arr( 1 ), angle( z_arr( 3, 1 ) ), 300, 'o', 'filled', 'markerfacecolor', 'w', 'markeredgecolor', c_arr{ 3 }, 'linewidth', 5 );


% Set up the video writer
outputVideo = VideoWriter( 'videos/synchronization_imag/single_limit_cycle.mp4', 'MPEG-4' );
outputVideo.FrameRate = 60;
open( outputVideo );

% Time per frame
numFrames    = round( T*outputVideo.FrameRate);
timePerFrame = T / numFrames;

for frameIdx = 1:numFrames

    % Current time for this frame
    currentTime = (frameIdx - 1) * timePerFrame;
    
    % Find the closest time in your timeArray (or interpolate if needed)
    [~,  idx ] = min( abs( t_arr - currentTime ) );

    set( p1, 'xdata',  real( z_arr( 1, idx ) ), 'ydata', imag( z_arr( 1, idx ) ) )
    set( p2, 'xdata',  real( z_arr( 2, idx ) ), 'ydata', imag( z_arr( 2, idx ) ) )
    set( p3, 'xdata',  real( z_arr( 3, idx ) ), 'ydata', imag( z_arr( 3, idx ) ) )

    set( p4, 'xdata',  t_arr( idx ), 'ydata', angle( z_arr( 1, idx ) ) )
    set( p5, 'xdata',  t_arr( idx ), 'ydata', angle( z_arr( 2, idx ) ) )
    set( p6, 'xdata',  t_arr( idx ), 'ydata', angle( z_arr( 3, idx ) ) )


    % Capture frame
    frame = getframe( f );
    
    % Write frame to video
    writeVideo( outputVideo, frame );
end

% Close the video file
close( outputVideo );

% Optionally, close the figure
close( f );


%% (1D) Andronov Hopf Oscillators, Single Limit Cycle, Revisit Pham

close all;

M1 = zeros( 3, 3 );

% Now, define three complex numbers for the three oscillators 
% Parameters of the Andronov-Hopf Oscillator 
r = 1.0; w = 2*pi;

% Using complex number
r0 = 0.8;
z1 = r0*exp(1i*(pi/2+0.1*pi)); 
z2 = r0*exp(1i*(pi/2-0.1*pi)); 
z3 = r0*exp(1i*pi/2); 

% First, run the simulation for both oscillators
dt = 1e-3;
T  = 8.0;
N  = round( T/dt );

% Time array 
t_arr = (0:N)*dt;
z_arr = zeros( 3, N+1 );
z_arr( :, 1 ) = [ z1; z2; z3 ];
M = M1;

k = 3;
M2 = k * [ -1, exp( 1i*(2*pi/3) ),  0;
           0, -1, exp( 1i*(2*pi/3) );
           exp( 1i*(2*pi/3) ), 0, -1 ];

for i = 2 : N+1
    dz1 = ( (r^2-abs(z1)^2 ) + 1i*w ) * z1 + M( 1, : )*[z1;z2;z3];
    dz2 = ( (r^2-abs(z2)^2 ) + 1i*w ) * z2 + M( 2, : )*[z1;z2;z3];
    dz3 = ( (r^2-abs(z3)^2 ) + 1i*w ) * z3 + M( 3, : )*[z1;z2;z3];

    if i == round( N/4 )
        M=M2;
    end

    z1 = z1 + dz1*dt;
    z2 = z2 + dz2*dt;
    z3 = z3 + dz3*dt;

    z_arr( 1, i ) = z1;
    z_arr( 2, i ) = z2;
    z_arr( 3, i ) = z3;
end

% Do the subplot
f = figure( ); 
a1 = subplot( 2, 3, 1 ); set( a1, 'parent', f ); hold on; axis equal;
a2 = subplot( 2, 3, 2 ); set( a2, 'parent', f ); hold on; axis equal;
a3 = subplot( 2, 3, 3 ); set( a3, 'parent', f ); hold on; axis equal;
lw = 2;
set( a1, 'xlim', [-lw, lw], 'ylim', [-lw, lw] )
set( a2, 'xlim', [-lw, lw], 'ylim', [-lw, lw] )
set( a3, 'xlim', [-lw, lw], 'ylim', [-lw, lw] )

plot( a1, real( z_arr( 1, : ) ), imag( z_arr( 1, : ) ), 'linewidth', 5, 'color', c_arr{ 1 } );
plot( a2, real( z_arr( 2, : ) ), imag( z_arr( 2, : ) ), 'linewidth', 5, 'color', c_arr{ 2 } );
plot( a3, real( z_arr( 3, : ) ), imag( z_arr( 3, : ) ), 'linewidth', 5, 'color', c_arr{ 3 } );

p1 = scatter( a1, real( z_arr( 1, 1 ) ), imag( z_arr( 1, 1 ) ), 300, 'o', 'filled', 'markerfacecolor', 'w', 'markeredgecolor', c_arr{ 1 }, 'linewidth', 5 );
p2 = scatter( a2, real( z_arr( 2, 1 ) ), imag( z_arr( 2, 1 ) ), 300, 'o', 'filled', 'markerfacecolor', 'w', 'markeredgecolor', c_arr{ 2 }, 'linewidth', 5 );
p3 = scatter( a3, real( z_arr( 3, 1 ) ), imag( z_arr( 3, 1 ) ), 300, 'o', 'filled', 'markerfacecolor', 'w', 'markeredgecolor', c_arr{ 3 }, 'linewidth', 5 );


a4 = subplot( 2, 3, 4:6 ); set( a4, 'parent', f ); hold on; 
set( a4, 'xlim', [0, max(t_arr)] )
plot( a4, t_arr, angle( z_arr( 1, : ) ), 'linewidth', 5, 'color', c_arr{ 1 } );
plot( a4, t_arr, angle( z_arr( 2, : ) ), 'linewidth', 5, 'color', c_arr{ 2 } );
plot( a4, t_arr, angle( z_arr( 3, : ) ), 'linewidth', 5, 'color', c_arr{ 3 } );
p4 = scatter( a4, t_arr( 1 ), angle( z_arr( 1, 1 ) ), 300, 'o', 'filled', 'markerfacecolor', 'w', 'markeredgecolor', c_arr{ 1 }, 'linewidth', 5 );
p5 = scatter( a4, t_arr( 1 ), angle( z_arr( 2, 1 ) ), 300, 'o', 'filled', 'markerfacecolor', 'w', 'markeredgecolor', c_arr{ 2 }, 'linewidth', 5 );
p6 = scatter( a4, t_arr( 1 ), angle( z_arr( 3, 1 ) ), 300, 'o', 'filled', 'markerfacecolor', 'w', 'markeredgecolor', c_arr{ 3 }, 'linewidth', 5 );


% Set up the video writer
outputVideo = VideoWriter( 'videos/synchronization_imag/two_limit_cycles.mp4', 'MPEG-4' );
outputVideo.FrameRate = 60;
open( outputVideo );

% Time per frame
numFrames    = round( T*outputVideo.FrameRate);
timePerFrame = T / numFrames;

for frameIdx = 1:numFrames

    % Current time for this frame
    currentTime = (frameIdx - 1) * timePerFrame;
    
    % Find the closest time in your timeArray (or interpolate if needed)
    [~,  idx ] = min( abs( t_arr - currentTime ) );

    set( p1, 'xdata',  real( z_arr( 1, idx ) ), 'ydata', imag( z_arr( 1, idx ) ) )
    set( p2, 'xdata',  real( z_arr( 2, idx ) ), 'ydata', imag( z_arr( 2, idx ) ) )
    set( p3, 'xdata',  real( z_arr( 3, idx ) ), 'ydata', imag( z_arr( 3, idx ) ) )

    set( p4, 'xdata',  t_arr( idx ), 'ydata', angle( z_arr( 1, idx ) ) )
    set( p5, 'xdata',  t_arr( idx ), 'ydata', angle( z_arr( 2, idx ) ) )
    set( p6, 'xdata',  t_arr( idx ), 'ydata', angle( z_arr( 3, idx ) ) )


    % Capture frame
    frame = getframe( f );
    
    % Write frame to video
    writeVideo( outputVideo, frame );
end

% Close the video file
close( outputVideo );

% Optionally, close the figure
close( f );

%% (1E) Andronov Hopf Oscillators, Parallel Combination of Synchronization

% Now, define three complex numbers for the three oscillators 
% Parameters of the Andronov-Hopf Oscillator 
r = 1.0; w = 2*pi;

% Using complex number
r0 = 0.8;
z1 = r0*exp(1i*(pi/2+2/3*pi)); 
z2 = r0*exp(1i*(pi/2-3/4*pi)); 
z3 = r0*exp(1i*pi/2); 

% First, run the simulation for both oscillators
dt = 1e-3;
T  = 8.0;
N  = round( T/dt );

% Time array 
t_arr = (0:N)*dt;
z_arr = zeros( 3, N+1 );
z_arr( :, 1 ) = [ z1; z2; z3 ];

% This is a zero matrix without any synchronization
M1 = zeros( 3, 3 );

% Two matrices that satisfies the row-sum-0 condition
M2 = [ -2,  1,  1; 
        1, -2,  1;
        1 , 1, -2];


M3 = [  -2,  1+1i,  1-1i; 
         1-1i, -2,  1+1i;
         1+1i,  1-1i, -2];

M = M1;

for i = 2 : N+1
    dz1 = ( (r^2-abs(z1)^2 ) + 1i*w ) * z1 + M( 1, : )*[z1;z2;z3];
    dz2 = ( (r^2-abs(z2)^2 ) + 1i*w ) * z2 + M( 2, : )*[z1;z2;z3];
    dz3 = ( (r^2-abs(z3)^2 ) + 1i*w ) * z3 + M( 3, : )*[z1;z2;z3];

    if i == round( N/4 )
        % M=M2;
        % M=M3;
        M = 3*M2 + 2*M3;
    end

    z1 = z1 + dz1*dt;
    z2 = z2 + dz2*dt;
    z3 = z3 + dz3*dt;

    z_arr( 1, i ) = z1;
    z_arr( 2, i ) = z2;
    z_arr( 3, i ) = z3;
end

% Do the subplot
f = figure( ); 
a1 = subplot( 2, 3, 1 ); set( a1, 'parent', f ); hold on; axis equal;
a2 = subplot( 2, 3, 2 ); set( a2, 'parent', f ); hold on; axis equal;
a3 = subplot( 2, 3, 3 ); set( a3, 'parent', f ); hold on; axis equal;
lw = 2;
set( a1, 'xlim', [-lw, lw], 'ylim', [-lw, lw] )
set( a2, 'xlim', [-lw, lw], 'ylim', [-lw, lw] )
set( a3, 'xlim', [-lw, lw], 'ylim', [-lw, lw] )

plot( a1, real( z_arr( 1, : ) ), imag( z_arr( 1, : ) ), 'linewidth', 5, 'color', c_arr{ 1 } );
plot( a2, real( z_arr( 2, : ) ), imag( z_arr( 2, : ) ), 'linewidth', 5, 'color', c_arr{ 2 } );
plot( a3, real( z_arr( 3, : ) ), imag( z_arr( 3, : ) ), 'linewidth', 5, 'color', c_arr{ 3 } );

p1 = scatter( a1, real( z_arr( 1, 1 ) ), imag( z_arr( 1, 1 ) ), 300, 'o', 'filled', 'markerfacecolor', 'w', 'markeredgecolor', c_arr{ 1 }, 'linewidth', 5 );
p2 = scatter( a2, real( z_arr( 2, 1 ) ), imag( z_arr( 2, 1 ) ), 300, 'o', 'filled', 'markerfacecolor', 'w', 'markeredgecolor', c_arr{ 2 }, 'linewidth', 5 );
p3 = scatter( a3, real( z_arr( 3, 1 ) ), imag( z_arr( 3, 1 ) ), 300, 'o', 'filled', 'markerfacecolor', 'w', 'markeredgecolor', c_arr{ 3 }, 'linewidth', 5 );


a4 = subplot( 2, 3, 4:6 ); set( a4, 'parent', f ); hold on; 
set( a4, 'xlim', [0, max(t_arr)] )
plot( a4, t_arr, angle( z_arr( 1, : ) ), 'linewidth', 5, 'color', c_arr{ 1 } );
plot( a4, t_arr, angle( z_arr( 2, : ) ), 'linewidth', 5, 'color', c_arr{ 2 } );
plot( a4, t_arr, angle( z_arr( 3, : ) ), 'linewidth', 5, 'color', c_arr{ 3 } );
p4 = scatter( a4, t_arr( 1 ), angle( z_arr( 1, 1 ) ), 300, 'o', 'filled', 'markerfacecolor', 'w', 'markeredgecolor', c_arr{ 1 }, 'linewidth', 5 );
p5 = scatter( a4, t_arr( 1 ), angle( z_arr( 2, 1 ) ), 300, 'o', 'filled', 'markerfacecolor', 'w', 'markeredgecolor', c_arr{ 2 }, 'linewidth', 5 );
p6 = scatter( a4, t_arr( 1 ), angle( z_arr( 3, 1 ) ), 300, 'o', 'filled', 'markerfacecolor', 'w', 'markeredgecolor', c_arr{ 3 }, 'linewidth', 5 );


% Set up the video writer
outputVideo = VideoWriter( 'videos/synchronization_imag/two_limit_cycles.mp4', 'MPEG-4' );
outputVideo.FrameRate = 60;
open( outputVideo );

% Time per frame
numFrames    = round( T*outputVideo.FrameRate);
timePerFrame = T / numFrames;

for frameIdx = 1:numFrames

    % Current time for this frame
    currentTime = (frameIdx - 1) * timePerFrame;
    
    % Find the closest time in your timeArray (or interpolate if needed)
    [~,  idx ] = min( abs( t_arr - currentTime ) );

    set( p1, 'xdata',  real( z_arr( 1, idx ) ), 'ydata', imag( z_arr( 1, idx ) ) )
    set( p2, 'xdata',  real( z_arr( 2, idx ) ), 'ydata', imag( z_arr( 2, idx ) ) )
    set( p3, 'xdata',  real( z_arr( 3, idx ) ), 'ydata', imag( z_arr( 3, idx ) ) )

    set( p4, 'xdata',  t_arr( idx ), 'ydata', angle( z_arr( 1, idx ) ) )
    set( p5, 'xdata',  t_arr( idx ), 'ydata', angle( z_arr( 2, idx ) ) )
    set( p6, 'xdata',  t_arr( idx ), 'ydata', angle( z_arr( 3, idx ) ) )


    % Capture frame
    frame = getframe( f );
    
    % Write frame to video
    writeVideo( outputVideo, frame );
end

% Close the video file
close( outputVideo );

% Optionally, close the figure
close( f );



%% (1E) Andronov Hopf Oscillators, Parallel Combination of Synchronization, Phase Offset

% Now, define three complex numbers for the three oscillators 
% Parameters of the Andronov-Hopf Oscillator 
r = 1.0; w = 2*pi;

% Using complex number
r0 = 0.8;
z1 = r0*exp(1i*(pi/2+1/10*pi)); 
z2 = r0*exp(1i*(pi/2-1/10*pi)); 
z3 = r0*exp(1i*pi/2); 

% First, run the simulation for both oscillators
dt = 1e-3;
T  = 8.0;
N  = round( T/dt );

% Time array 
t_arr = (0:N)*dt;
z_arr = zeros( 3, N+1 );
z_arr( :, 1 ) = [ z1; z2; z3 ];

% This is a zero matrix without any synchronization
M1 = zeros( 3, 3 );

% Two matrices that satisfies the row-sum-0 condition
k = 3;
M2 = k * [ -1,  exp( 1i*2*pi/3 ),  0; 
            0, -1,  exp( 1i*2*pi/3 );
            exp( 1i*2*pi/3 ) , 0, -1];


M3 = k*[  -1, exp( 1i*2*pi/3 )*0.5,  exp( 1i*4*pi/3 )*0.5; 
            exp( 1i*4*pi/3 )*0.5, -1,  exp( 1i*2*pi/3 )*0.5;
            exp( 1i*2*pi/3 )*0.5,  exp( 1i*4*pi/3 )*0.5, -1];

M = M1;

for i = 2 : N+1
    dz1 = ( (r^2-abs(z1)^2 ) + 1i*w ) * z1 + M( 1, : )*[z1;z2;z3];
    dz2 = ( (r^2-abs(z2)^2 ) + 1i*w ) * z2 + M( 2, : )*[z1;z2;z3];
    dz3 = ( (r^2-abs(z3)^2 ) + 1i*w ) * z3 + M( 3, : )*[z1;z2;z3];

    if i == round( N/4 )
        M = M3;
        % M = 0.4*M2 + 1*M3;
    end

    z1 = z1 + dz1*dt;
    z2 = z2 + dz2*dt;
    z3 = z3 + dz3*dt;

    z_arr( 1, i ) = z1;
    z_arr( 2, i ) = z2;
    z_arr( 3, i ) = z3;
end

% Do the subplot
f = figure( ); 
a1 = subplot( 2, 3, 1 ); set( a1, 'parent', f ); hold on; axis equal;
a2 = subplot( 2, 3, 2 ); set( a2, 'parent', f ); hold on; axis equal;
a3 = subplot( 2, 3, 3 ); set( a3, 'parent', f ); hold on; axis equal;
lw = 2;
set( a1, 'xlim', [-lw, lw], 'ylim', [-lw, lw] )
set( a2, 'xlim', [-lw, lw], 'ylim', [-lw, lw] )
set( a3, 'xlim', [-lw, lw], 'ylim', [-lw, lw] )

plot( a1, real( z_arr( 1, : ) ), imag( z_arr( 1, : ) ), 'linewidth', 5, 'color', c_arr{ 1 } );
plot( a2, real( z_arr( 2, : ) ), imag( z_arr( 2, : ) ), 'linewidth', 5, 'color', c_arr{ 2 } );
plot( a3, real( z_arr( 3, : ) ), imag( z_arr( 3, : ) ), 'linewidth', 5, 'color', c_arr{ 3 } );

p1 = scatter( a1, real( z_arr( 1, 1 ) ), imag( z_arr( 1, 1 ) ), 300, 'o', 'filled', 'markerfacecolor', 'w', 'markeredgecolor', c_arr{ 1 }, 'linewidth', 5 );
p2 = scatter( a2, real( z_arr( 2, 1 ) ), imag( z_arr( 2, 1 ) ), 300, 'o', 'filled', 'markerfacecolor', 'w', 'markeredgecolor', c_arr{ 2 }, 'linewidth', 5 );
p3 = scatter( a3, real( z_arr( 3, 1 ) ), imag( z_arr( 3, 1 ) ), 300, 'o', 'filled', 'markerfacecolor', 'w', 'markeredgecolor', c_arr{ 3 }, 'linewidth', 5 );


a4 = subplot( 2, 3, 4:6 ); set( a4, 'parent', f ); hold on; 
set( a4, 'xlim', [0, max(t_arr)] )
plot( a4, t_arr, angle( z_arr( 1, : ) ), 'linewidth', 5, 'color', c_arr{ 1 } );
plot( a4, t_arr, angle( z_arr( 2, : ) ), 'linewidth', 5, 'color', c_arr{ 2 } );
plot( a4, t_arr, angle( z_arr( 3, : ) ), 'linewidth', 5, 'color', c_arr{ 3 } );
p4 = scatter( a4, t_arr( 1 ), angle( z_arr( 1, 1 ) ), 300, 'o', 'filled', 'markerfacecolor', 'w', 'markeredgecolor', c_arr{ 1 }, 'linewidth', 5 );
p5 = scatter( a4, t_arr( 1 ), angle( z_arr( 2, 1 ) ), 300, 'o', 'filled', 'markerfacecolor', 'w', 'markeredgecolor', c_arr{ 2 }, 'linewidth', 5 );
p6 = scatter( a4, t_arr( 1 ), angle( z_arr( 3, 1 ) ), 300, 'o', 'filled', 'markerfacecolor', 'w', 'markeredgecolor', c_arr{ 3 }, 'linewidth', 5 );


% Set up the video writer
outputVideo = VideoWriter( 'videos/synchronization_imag/two_limit_cycles.mp4', 'MPEG-4' );
outputVideo.FrameRate = 60;
open( outputVideo );

% Time per frame
numFrames    = round( T*outputVideo.FrameRate);
timePerFrame = T / numFrames;

for frameIdx = 1:numFrames

    % Current time for this frame
    currentTime = (frameIdx - 1) * timePerFrame;
    
    % Find the closest time in your timeArray (or interpolate if needed)
    [~,  idx ] = min( abs( t_arr - currentTime ) );

    set( p1, 'xdata',  real( z_arr( 1, idx ) ), 'ydata', imag( z_arr( 1, idx ) ) )
    set( p2, 'xdata',  real( z_arr( 2, idx ) ), 'ydata', imag( z_arr( 2, idx ) ) )
    set( p3, 'xdata',  real( z_arr( 3, idx ) ), 'ydata', imag( z_arr( 3, idx ) ) )

    set( p4, 'xdata',  t_arr( idx ), 'ydata', angle( z_arr( 1, idx ) ) )
    set( p5, 'xdata',  t_arr( idx ), 'ydata', angle( z_arr( 2, idx ) ) )
    set( p6, 'xdata',  t_arr( idx ), 'ydata', angle( z_arr( 3, idx ) ) )


    % Capture frame
    frame = getframe( f );
    
    % Write frame to video
    writeVideo( outputVideo, frame );
end

% Close the video file
close( outputVideo );

% Optionally, close the figure
close( f );
