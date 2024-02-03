% [Title]     Study of Synchronization of Multiple Limit Cycles
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


%% (1A) Plot the Andronov Hopf Oscillator 
close all;

is_mode = 2;

Nosc = 4;
assert( ismember(    Nosc, 1:Nc ) );
assert( ismember( is_mode, 1:2  ) );

osc_arr   = cell( 1, Nosc );
x_arr     = cell( 1, Nosc );
phase_arr = cell( 1, Nosc );

% Radius and angular velcoity
r = 1.0;
w = pi;

for i = 1 : Nosc
    osc_arr{ i } = AndronovHopfOSC( w, r );
end

% First, run the simulation for both oscillators
dt = 1e-3;
T  = 8.0;
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
        x_new = osc.step( dt, x( :, i-1), zeros( 2, 1 ) );
        
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
    
    if is_mode == 1
        a2 = subplot( 2, Nosc, Nosc+i   );
        hold( a2, 'on' ); 

        plot( a2, t_arr, x( 1, : ), 'linewidth', 3, 'color', 'r' )
        plot( a2, t_arr, x( 2, : ), 'linewidth', 3, 'color', 'g' )

        % Define an array for 2 elements on the second part
        p_markers{ 2, i } = zeros( 1, 2 );
        p_markers{ 2, i }( 1 ) = scatter( a2, t_arr( 1 ), x( 1, 1 ), m_size, 'markeredgecolor', 'r', 'markerfacecolor', 'w', 'linewidth', 4  );
        p_markers{ 2, i }( 2 ) = scatter( a2, t_arr( 1 ), x( 2, 1 ), m_size, 'markeredgecolor', 'g', 'markerfacecolor', 'w', 'linewidth', 4  );
        set( a2, 'ylim', [-lw,lw], 'xlim', [0, max( t_arr ) ] )

    elseif is_mode == 2
        a2 = subplot( 2, Nosc, [Nosc+1, 2*Nosc]   );
        hold( a2, 'on' ); 
        plot( a2, t_arr, phase_arr{ i }, 'color', c_arr{ i }  )
        p_markers{ 2, i } = scatter( a2, t_arr( 1 ), phase_arr{ i }( 1 ), m_size, 'markeredgecolor', c_arr{ i }, 'markerfacecolor', 'w', 'linewidth', 4  );
        set( a2, 'ylim', [-0.1,2.1]*pi, 'xlim', [0, max( t_arr ) ] )
    end

    p_markers{ 1, i } = scatter( a1, x( 1, 1 ), x( 2, 1 ), m_size, 'markerfacecolor', 'w', 'markeredgecolor', c_arr{ i }, 'linewidth', 4  );

end

% Set up the video writer
outputVideo = VideoWriter( ['videos/synchronization/simple_osc', num2str( is_mode ), '.mp4'], 'MPEG-4' );
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
         set( p_markers{ 1, i }     , 'xdata',  x( 1, idx ), 'ydata', x( 2, idx ) )
         if is_mode == 1
             set( p_markers{ 2, i }( 1 ), 'xdata', t_arr( idx ), 'ydata', x( 1, idx ) )
             set( p_markers{ 2, i }( 2 ), 'xdata', t_arr( idx ), 'ydata', x( 2, idx ) )

         elseif is_mode == 2
             set( p_markers{ 2, i }, 'xdata', t_arr( idx ), 'ydata', phase_arr{ i }( idx ) )
         else

         end
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

%% (1B) Synchronization, Excitator Only
close all;
% The k should be high enough for synchronization
% Section 5.3 of Pham, Quang-Cuong, and Jean-Jacques Slotine. 
% "Stable concurrent synchronization in dynamic system networks." (2007).

Nosc = 6;
osc_arr   = cell( 1, Nosc );
x_arr     = cell( 1, Nosc );
phase_arr = cell( 1, Nosc );


% Radius and angular velcoity
r = 1.0;
w = pi;

for i = 1 : Nosc
    osc_arr{ i } = AndronovHopfOSC( w, r );
end


% First, run the simulation for both oscillators
dt = 1e-4;
T  = 6.0;
N  = round( T/dt );

% Time array 
t_arr = (0:N)*dt;

% Initialize the array of x
for i = 1 : Nosc
    x_arr{ i } = zeros( 2, N+1 );
    x_arr{ i }( :, 1 ) = 2*rand( 2, 1 )-2;
end

% The angle array 
ang_arr = 2*pi/Nosc * ones( 1, Nosc );

% Gain
k = 0;

for i = 2 : N+1

    if i == round( N/3 )
        k = 1e-3;
    end

    % Iterate through the oscillator
    for j = 1 : Nosc
        osc = osc_arr{ j };
        x   = x_arr{ j };
        ang = ang_arr( j );

        if j == Nosc
            tmp = 1;
        else
            tmp = j+1;
        end
        input = k* [ cos( ang ), -sin( ang );
                     sin( ang ),  cos( ang ) ] * x_arr{ tmp }( :, i-1 );

        x_new = osc.step( dt, x( :, i-1), input );
        
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
f = figure( ); 

% Draw the background image
% Set the markers of each image 
p_markers = cell( 2, Nosc );

m_size = 200;

for i = 1 : Nosc
    x = x_arr{ i };
    a1 = subplot( 2, Nosc, i ); 
    axis equal
    a2 = subplot( 2, Nosc, [Nosc+1:2*Nosc]   );
    hold( a1, 'on' ); hold( a2, 'on' ); 
    set( a2, 'ylim', [-0.1,2.1]*pi, 'xlim', [0, max( t_arr ) ] )
   
    plot( a1, x( 1, : ), x( 2, : ), 'linewidth', 3, 'color', c_arr{ i }, 'linestyle', '-' )
    plot( a2, t_arr, phase_arr{ i }, 'linewidth', 3, 'color', c_arr{ i }, 'linestyle', '-' )    
    
    p_markers{ 1, i } = scatter( a1, x( 1, 1 ), x( 2, 1 ), m_size, 'markerfacecolor', 'w', 'markeredgecolor', c_arr{ i }, 'linewidth', 4  );
    p_markers{ 2, i } = scatter( a2, t_arr( 1 ), phase_arr{ i }( 1 ), m_size, 'markerfacecolor', 'w', 'markeredgecolor', c_arr{ i }, 'linewidth', 4  );
end


% Set up the video writer
outputVideo = VideoWriter( 'videos/synchronization/excitator.mp4', 'MPEG-4' );
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
         set( p_markers{ 2, i }, 'xdata',  t_arr( idx ), 'ydata', phase_arr{ i }( idx ) )
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


%% (1C) Synchronization, Diffusive Model
close all;
% The k should be high enough for synchronization
% Section 5.3 of Pham, Quang-Cuong, and Jean-Jacques Slotine. 
% "Stable concurrent synchronization in dynamic system networks." (2007).

Nosc = 4;
ang_del = 2*pi/Nosc;
osc_arr = cell( 1, Nosc );
x_arr   = cell( 1, Nosc );

% Radius and angular velcoity
r = 1.0;
w = pi;

for i = 1 : Nosc
    osc_arr{ i } = AndronovHopfOSC( w, r );
end

% First, run the simulation for both oscillators
dt = 1e-4;
T  = 8.0;
N  = round( T/dt );

% Time array 
t_arr = (0:N)*dt;

% Initialize the array of x
for i = 1 : Nosc
    x_arr{ i } = zeros( 2, N+1 );
    x_arr{ i }( :, 1 ) = 2*rand( 2, 1 )-2;
end

% Gain
k = 0;

for i = 2 : N+1

    if i == round( 0.25*N )
        k = 1e-3;
    end

    if i == round( 0.7*N )
        for j = 1 : Nosc
            osc_arr{ j }.w = 2*pi;
        end
    end

    % Iterate through the oscillator
    for j = 1 : Nosc
        osc = osc_arr{ j };
        x   = x_arr{ j };

        if j == Nosc
            tmp = 1;
        else
            tmp = j+1;
        end
        input = k* [ cos( ang_del ), -sin( ang_del );
                     sin( ang_del ),  cos( ang_del ) ] * x_arr{ tmp }( :, i-1 ) - k * x_arr{ j }( :, i-1 );

        x_new = osc.step( dt, x( :, i-1), input );
        
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
f = figure( ); 

% Draw the background image
% Set the markers of each image 
p_markers = cell( 2, Nosc );

m_size = 200;
lw = 2*r

for i = 1 : Nosc
    x = x_arr{ i };
    a1 = subplot( 2, Nosc, i ); 
    axis equal
    a2 = subplot( 2, Nosc, [Nosc+1:2*Nosc]   );
    hold( a1, 'on' ); hold( a2, 'on' ); 
    set( a1, 'xlim', [-lw, lw], 'ylim', [-lw, lw] )
    set( a2, 'ylim', [-0.1,2.1]*pi, 'xlim', [0, max( t_arr ) ] )
   
    plot( a1, x( 1, : ), x( 2, : ), 'linewidth', 3, 'color', c_arr{ i }, 'linestyle', '-' )
    plot( a2, t_arr, phase_arr{ i }, 'linewidth', 3, 'color', c_arr{ i }, 'linestyle', '-' )    
    
    p_markers{ 1, i } = scatter( a1, x( 1, 1 ), x( 2, 1 ), m_size, 'markerfacecolor', 'w', 'markeredgecolor', c_arr{ i }, 'linewidth', 4  );
    p_markers{ 2, i } = scatter( a2, t_arr( 1 ), phase_arr{ i }( 1 ), m_size, 'markerfacecolor', 'w', 'markeredgecolor', c_arr{ i }, 'linewidth', 4  );
end


% Set up the video writer
outputVideo = VideoWriter( 'videos/synchronization/excitator.mp4', 'MPEG-4' );
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
         set( p_markers{ 2, i }, 'xdata',  t_arr( idx ), 'ydata', phase_arr{ i }( idx ) )
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

