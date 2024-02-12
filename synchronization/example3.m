% [Title]     Study of Synchronization, Fithugh-Nagumo
% [Author]    Moses Chong-ook Nah
% [Email]     mosesnah@mit.edu
% [Update]    At 2024.02.11

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

%% (1A) Simple Oscillator


Nosc = 1;
assert( ismember(    Nosc, 1:Nc ) );

osc_arr   = cell( 1, Nosc );
x_arr     = cell( 1, Nosc );
phase_arr = cell( 1, Nosc );

a = 0.7;
b = 0.8;
c = 8.0;

for i = 1 : Nosc
    osc_arr{ i } = FitzHughNagumoOSC( a, b, c );
end

% First, run the simulation for both oscillators
dt = 1e-3;
T  = 50.0;
N  = round( T/dt );

% Time array 
t_arr = (0:N)*dt;

% Initialize the array of x
for i = 1 : Nosc
    x_arr{ i } = zeros( 2, N+1 );
    x_arr{ i }( :, 1 ) = 2*rand( 2, 1 )-2;
end

for i = 2 : N+1

    % Iterate through the oscillator
    for j = 1 : Nosc
        osc = osc_arr{ j };
        x   = x_arr{ j };
        x_new = osc.step( dt, x( :, i-1), -1.4, zeros( 2, 1 ) );
        
        x_arr{ j }( :, i ) = x_new;
    end
    
end

% The arc tangent
for i = 1 : Nosc
    x = x_arr{ i };
    x1 = x( 1, : );
    x2 = x( 2, : );
    phase_arr{ i } = mod(atan2( x2,x1 ),2*pi);
end

% Draw the video 
lw= 2.0;
f = figure( ); 

% Draw the background image
% Set the markers of each image 
p_markers = cell( 2, Nosc );

m_size = 200;

for i = 1 : Nosc
    x = x_arr{ i };
    a1 = subplot( 2, Nosc, i ); 
    axis equal
    hold( a1, 'on' )
    set( a1, 'xlim', [-lw,lw], 'ylim', [-lw,lw] )
    
    plot( a1, x( 1, : ), x( 2, : ), 'linewidth', 3, 'color', c_arr{ i } )
    

    a2 = subplot( 2, Nosc, [Nosc+1, 2*Nosc]   );
    hold( a2, 'on' ); 
    plot( a2, t_arr, phase_arr{ i }, 'color', c_arr{ i }  )
    p_markers{ 2, i } = scatter( a2, t_arr( 1 ), phase_arr{ i }( 1 ), m_size, 'markeredgecolor', c_arr{ i }, 'markerfacecolor', 'w', 'linewidth', 4  );
    set( a2, 'ylim', [-0.1,2.1]*pi, 'xlim', [0, max( t_arr ) ] )

    p_markers{ 1, i } = scatter( a1, x( 1, 1 ), x( 2, 1 ), m_size, 'markerfacecolor', 'w', 'markeredgecolor', c_arr{ i }, 'linewidth', 4  );

end

% Set up the video writer
outputVideo = VideoWriter( 'videos/synchronization/simple_osc.mp4', 'MPEG-4' );
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
   
    for i = 1 : Nosc
         x = x_arr{ i };
         set( p_markers{ 1, i }, 'xdata',  x( 1, idx ), 'ydata', x( 2, idx ) )
         set( p_markers{ 2, i }, 'xdata', t_arr( idx ), 'ydata', phase_arr{ i }( idx ) )
    end

    
    % Capture frame
    frame = getframe( f );
    
    % Write frame to video
    writeVideo( outputVideo, frame );
end

% Close the video file
close( outputVideo );

% Optionally, close the figure
close( f );
