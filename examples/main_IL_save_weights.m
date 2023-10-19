%% ==================================================================
%% [Title] Imitation Learning and Saved the Weight Array
% Author: Moses Chong-ook Nah
%  Email: mosesnah@mit.edu
%   Date: 2023.10.18
%% ==================================================================

%% [0A] Initialization
clear; close all; clc;
fig_config( 'fontSize', 20, 'markerSize', 10 )

%% [1-] Imitation Learning of a Given Trajectory
%% ---- [1A] Define Trajectory

is_check = true; 

% Define the Analytical Trajectory 2D that we aim to learn
syms t_sym

% Name of the Trajectory which we aim to learn
name_traj = 'cosine';

% Define the Trajectory, which is a 3D example.
p0i = [ 0, 0, 0 ];
D   = 3.0;

% Position
switch name_traj

    case 'min_jerk'
        p0f = [ 1, 1, 1 ];      
        tn  = t_sym/D;

        px = p0i( 1 ) + ( p0f( 1 ) - p0i( 1 ) ) * ( 10 * tn^3 - 15 * tn^4 + 6 * tn^5 );
        py = p0i( 2 ) + ( p0f( 2 ) - p0i( 2 ) ) * ( 10 * tn^3 - 15 * tn^4 + 6 * tn^5 );
        pz = p0i( 3 ) + ( p0f( 3 ) - p0i( 3 ) ) * ( 10 * tn^3 - 15 * tn^4 + 6 * tn^5 );

    case 'cosine'
        l  = 1.5;
        D  = 5.0;
        v  = l/D;        
        tn = t_sym/D;        
        px = p0i( 1 );
        py = p0i( 2 ) + l * ( 10 * tn^3 - 15 * tn^4 + 6 * tn^5 );
        pz = p0i( 3 ) + 2*( cos( 2*2*pi/l * ( l * tn ) ) - 1 );
    
        % Add some rotation too
        p_tmp = [px; py; pz];
        p_tmp = rotx( 20 ) * p_tmp;
        p_tmp = roty( 30 ) * p_tmp;
        p_tmp = rotz( 20 ) * p_tmp;

        px = p_tmp( 1 );
        py = p_tmp( 2 );
        pz = p_tmp( 3 );

        p0f = double( subs( p_tmp, t_sym, D )' );

    case 'radial'
        r0i = 3.0;
        r0f = 1.0;
        D   = 5.0;

        tn = t_sym/D;         

        rt    = r0i + ( r0f - r0i ) * ( -2*tn^3 + 3*tn^2 );
        theta = 2 * pi * ( -tn^4 + 2*tn^2 );

        px = p0i( 1 );
        py = rt * cos( theta );
        pz = rt * sin( theta );

        % Add some rotation too
        p_tmp = [px; py; pz];
        p_tmp = rotx( 20 ) * p_tmp;
        p_tmp = roty( 30 ) * p_tmp;
        p_tmp = rotz( 20 ) * p_tmp;

        px = p_tmp( 1 );
        py = p_tmp( 2 );
        pz = p_tmp( 3 );

        p0i = double( subs( p_tmp, t_sym, 0 )' );
        p0f = double( subs( p_tmp, t_sym, D )' );

    otherwise

end
    

% Velocity
dpx = diff( px, t_sym );
dpy = diff( py, t_sym );
dpz = diff( pz, t_sym );

% Acceleration
ddpx = diff( dpx, t_sym );
ddpy = diff( dpy, t_sym );
ddpz = diff( dpz, t_sym );

% Saving these elements as symbolic array.
p_sym   = {   px,   py,   pz };
dp_sym  = {  dpx,  dpy,  dpz };
ddp_sym = { ddpx, ddpy, ddpz };

% Generating the actual data from the symbolic form 
% These data will be used for Imitation Learning
dt    = 0.01;
t_arr = 0:dt:D;
P     = length( t_arr );

  p_data = zeros( 3, P );
 dp_data = zeros( 3, P );
ddp_data = zeros( 3, P );

for i = 1 : P
    for j = 1 : 3
         p_data( j, i ) = double( subs(   p_sym{ j }, t_sym, t_arr( i ) ) );
        dp_data( j, i ) = double( subs(  dp_sym{ j }, t_sym, t_arr( i ) ) );
       ddp_data( j, i ) = double( subs( ddp_sym{ j }, t_sym, t_arr( i ) ) );
    end
end

if is_check 
    f = figure( );
    p1 = subplot( 3, 1, 1 );
    plot( p1, t_arr,   p_data );

    p2 = subplot( 3, 1, 2 );
    plot( p2, t_arr,  dp_data );
    
    p3 = subplot( 3, 1, 3 );
    plot( p3, t_arr,  ddp_data );
    
end

f = figure( ); a = axes( 'parent', f );
plot3( a, p_data( 1, : ), p_data( 2, : ), p_data( 3, : ))

%% ---- [1B] Learning the Weights for Imitation Learning

% For Imitation Learning, one should define the number of basis function
N  = 50;

% Parameters of the 3 DMPs
% We use the identical alpha_z, beta_z values
alpha_z = 10.0;
alpha_s = 1.0;
beta_z  = 1/4 * alpha_z;
tau     = D;
z0      = 0;

cs  = CanonicalSystem( 'discrete', tau, alpha_s );

tsx = TransformationSystem( alpha_z, beta_z, tau, p_data( 1, 1 ), dp_data( 1, 1 )*tau );
tsy = TransformationSystem( alpha_z, beta_z, tau, p_data( 2, 1 ), dp_data( 2, 1 )*tau );
tsz = TransformationSystem( alpha_z, beta_z, tau, p_data( 3, 1 ), dp_data( 3, 1 )*tau );
ts_arr = { tsx, tsy, tsz };

fs = NonlinearForcingTerm( cs, N );


%% ---- [1C] Learning Weights via Least-Square

w_arr = zeros( 3, N );

for i = 1 : 3

    if p0f( i ) - p0i( i ) == 0
        continue
    end

    phi_mat = zeros( P, N );    
    f_arr = ts_arr{ i }.get_desired( p_data( i, : ), dp_data( i, : ), ddp_data( i, : ), p0f( i ) );

    % Over the Time Array
    for j = 1 : P
        phi_sum = fs.calc_whole_at_t( t_arr( j ) );

        % Over the time basis function
        for k = 1 : N
            phi_mat( j, k ) = fs.calc_ith( t_arr( j ), k ) / phi_sum * ( p0f( i ) - p0i( i ) ) * cs.calc( t_arr( j ) );
        end
    end

    % Get w_arr with Least square solution
    ttmp = ( phi_mat' * phi_mat )^(-1) * phi_mat' * f_arr';
    w_arr( i, : ) = ttmp';

end

%% ---- [1D] Rollout Generating a Full trajectory with the Transformation System

% Resetting the Transformation System
for i = 1 : 3
    ts_arr{ i }.reset( )
end

% Initial Time and time step
t0i = 2.0;
dt  = 1e-3;
Nt  = 8000;

% The total time and its time array
T      = dt * Nt;
t_arr2 = dt * (0:(Nt-1));

% For plotting and saving the data
y_arr  = zeros( 3, Nt );
z_arr  = zeros( 3, Nt );
dy_arr = zeros( 3, Nt );
dz_arr = zeros( 3, Nt );

% Initial Condition
 y_arr( :, 1 )  =  p_data( :, 1 );
dy_arr( :, 1 )  = dp_data( :, 1 );
 z_arr( :, 1 )  = dp_data( :, 1 ) * tau;

% The analytical data 
  p_arr_test = zeros( 3, Nt );
 dp_arr_test = zeros( 3, Nt );
ddp_arr_test = zeros( 3, Nt );

t = 0;


for i = 0 : (Nt-1)
    
    for k = 1 : 3

        if sum( abs( w_arr( k, : ) ) ) == 0
            continue;
        end

        % Before conducting the movement
        % Maintaining that posture
        if t <= t0i
            y_arr( k, i + 1 ) =  y_arr( k, 1 );
           dy_arr( k, i + 1 ) = dy_arr( k, 1 );            
            z_arr( k, i + 1 ) =  z_arr( k, 1 );
    
              p_arr_test( k, i + 1 ) =   p_data( k, 1 );
             dp_arr_test( k, i + 1 ) =  dp_data( k, 1 );
            ddp_arr_test( k, i + 1 ) = ddp_data( k, 1 );
    
        % During the movement
        elseif t0i <= t
            
            if t<= t0i + D
    
                % taking off the initial time offset
                t_tmp = t - t0i;
    
                % Calculating the input from the weights
                % First, check if whole activation value is 0
                phi_sum = fs.calc_whole_at_t( t_tmp );
                
                f_input = 0;
    
                if phi_sum ~= 0
                    f_input = fs.calc_whole_weighted_at_t( t_tmp, w_arr( k, : ) )/phi_sum;
                    f_input = f_input*( p0f( k )- p0i( k ) )*cs.calc( t_tmp );
                end
                
                % Position Array
                  p_arr_test( k, i + 1 ) = double( subs(   p_sym{ k }, t_sym, t_tmp ) );
                 dp_arr_test( k, i + 1 ) = double( subs(  dp_sym{ k }, t_sym, t_tmp ) );
                ddp_arr_test( k, i + 1 ) = double( subs( ddp_sym{ k }, t_sym, t_tmp ) );       
    
            else
                f_input = 0; 
    
                  p_arr_test( k, i + 1 ) =   p_data( k, end );
                 dp_arr_test( k, i + 1 ) =    0;
                ddp_arr_test( k, i + 1 ) =    0;
            end
    
            [ y, z, dy, dz ] = ts_arr{ k }.step( p0f( k ), f_input, dt );
            y_arr(  k, i + 1 ) = y;
            z_arr(  k, i + 1 ) = z; 
            dy_arr( k, i + 1 ) = dy;
            dz_arr( k, i + 1 ) = dz;
    
        end

    end

    t = t + dt;
end

%% [1E] Plotting the Results
% Plot the results
f = figure( );
subplot( 3, 1, 1 )
hold on
plot( t_arr2, y_arr, 'linewidth', 3 )
plot( t_arr2, p_arr_test, 'linewidth', 3, 'linestyle', '--' )
set( gca, 'xlim', [0, T], 'fontsize', 30, 'xticklabel', {} )
ylabel( 'Pos. (-)' )

subplot( 3, 1, 2 )
hold on
plot( t_arr2, dy_arr, 'linewidth', 3 )
plot( t_arr2, dp_arr_test, 'linewidth', 3, 'linestyle','--' )
set( gca, 'xlim', [0, T], 'fontsize', 30, 'xticklabel', {} )
ylabel( 'Vel. (-)' )

subplot( 3, 1, 3 )
hold on
plot( t_arr2, dz_arr/tau, 'linewidth', 3 )
plot( t_arr2, ddp_arr_test, 'linewidth', 3, 'linestyle','--' )
set( gca, 'xlim', [0, T], 'fontsize', 30, 'xtick', 0:0.5:T)
xlabel( 'Time (sec)' )
ylabel( 'Acc. (-)' )

figure( )
plot3( p_arr_test( 1, : ), p_arr_test( 2, : ), p_arr_test( 3, : ), 'color', 'k' )

%% [1F] Saving the Data

% The Weighting Matrix and additional Data
data = struct;
data.p_sym   =   p_sym;
data.dp_sym  =  dp_sym;
data.ddp_sym = ddp_sym;

data.name = name_traj;
data.alpha_z = alpha_z;
data.beta_z  =  beta_z;
data.alpha_s = alpha_s;

data.weight  = w_arr;
data.p0i = p0i;
data.p0f = p0f;
data.tau = tau;

data.cs = cs;
data.fs = fs;
data.ts_arr = ts_arr;

save( ['learned_parameters/', name_traj ,'.mat'], 'data' );
