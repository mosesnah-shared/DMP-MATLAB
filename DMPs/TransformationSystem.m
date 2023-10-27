classdef TransformationSystem < handle
    
    properties
        
        % Parameters of the transformation systems
        alpha_z;
        beta_z;
        tau;

        % The canonical system 
        cs; 
 
    end
    
    methods
        function obj = TransformationSystem( alpha_z, beta_z, cs )
            % ===========================================================================
            % Descriptions
            % ------------
            %    Constructor of Transformation System 
            %    trans_sys = TransformationSystem( alpha_z, beta_z, cs )
            %       
            %    It is a second-order linear system with an input
            %       tau * dy = z
            %       tau * dz = alpha_z * [ beta_z * (g - y) - z ] + input
            %    
            %    dy, dz, z can be arrays
            %
            % Authors           
            % -------           
            %   Moses C. Nah    mosesnah@mit.edu
            % 
            % Parameters
            % ----------
            %   (1) alpha_z - A Positive Gain
            % 
            %   (2) beta_z  - A Positive Gain
            %
            %   (3) cs      - Canonical System, which provides 
            %                 the duration of the system
            %
            % ===========================================================================
            
            assert( alpha_z > 0 && beta_z > 0 );
           
            obj.alpha_z = alpha_z;
            obj.beta_z  =  beta_z;
            obj.cs      = cs;
            obj.tau     = cs.tau;
           
        end
        
        function [ y_new, z_new, dy, dz ] = step( obj, y_old, z_old, g, input, dt )
            % ===========================================================================            
            % Descriptions
            % ------------
            %    [ y_new, z_new, dy, dz ] = step( obj, y_old, z_old, g, input, dt  )
            %       
            %    A forward integration given the following differential equation
            %       tau * dy = z_old
            %       tau * dz = alpha_z * [ beta_z * ( g - y_old ) - z_old ] + input
            %    
            %    A simple single-step integration
            %       y_new = y_old + dt * dy   
            %       z_new = z_old + dt * dz   
            % 
            % Authors           
            % -------           
            %   Moses C. Nah    mosesnah@mit.edu
            % 
            % Parameters
            % ----------
            %   (1) y_old - y_old, the position
            %             
            %   (2) z_old - z_old, the generalized velocity
            %             
            %   (3) g     - Goal position, can be modified
            % 
            %   (4) input - force, could include coupling term
            %
            %   (5) dt    - time step
            %
            % Returns
            % -------
            %   (1) y_new - new position that is integrated
            % 
            %   (2) z_new - new generalized velocity that is integrated
            %
            %   (3) dy    - the velocity, z/tau
            %
            %   (4) dz    - the acceleration
            % ===========================================================================
            
            % y_old, z_old, g, input must all be column vectors 
            assert( iscolumn( y_old ) && iscolumn( z_old ) && iscolumn( g ) && iscolumn( input ) )
            
            % The length must also be identical 
            assert( length( y_old ) == length( z_old )  )
            assert( length( y_old ) == length( g     )  )
            assert( length( y_old ) == length( input )  )
            
            % Time-step must be strictly positive 
            assert( isscalar( dt ) && dt > 0 )
         
            % Time-derivative of y, z 
            dy = z_old/obj.tau;
            dz = ( obj.alpha_z * obj.beta_z * ( g - y_old ) - obj.alpha_z * z_old + input )/obj.tau;
            
            % Take out the y_new, z_new
            y_new = dy * dt + y_old;
            z_new = dz * dt + z_old;
            
        end
        
        function f_des = get_desired( obj, y_des, dy_des, ddy_des, g )
            % ===========================================================================            
            % Descriptions
            % ------------
            %    f_des = get_desired( obj, y_des, dy_des, ddy_des, g )
            %       
            %    Calculate an Array of f_des, from the y_desired trajectory.
            % 
            % Authors           
            % -------           
            %   Moses C. Nah    mosesnah@mit.edu
            % 
            % Parameters
            % ----------
            %   (1)   y_des -     Position of Desired Trajectory, nxP array
            % 
            %   (2)  dy_des -     Velocity of Desired Trajectory, nxP array
            %
            %   (3) ddy_des - Acceleration of Desired Trajectory, nxP array
            %
            %   (4)       g - Goal location
            % ===========================================================================
                
            % Check whether y_des, dy_des, ddy_des are all same size 
            assert( all( size( y_des ) == size(  dy_des ) ) );
            assert( all( size( y_des ) == size( ddy_des ) ) );
            assert( iscolumn( g ) )

            % Check whether the size of goal is consistent 
            [ n, ~ ] = size( y_des );
            assert( length( g ) == n );

            % Calculate the f_des array
            f_des = obj.tau^2 * ddy_des + obj.alpha_z * obj.tau * dy_des + obj.alpha_z * obj.beta_z * ( y_des - g );

        end

        function [ y_arr, z_arr, dy_arr ] = rollout( obj, y0, z0, g, input_arr, t0i, t_arr )
            % ===========================================================================            
            % Conduct integration 
            % ===========================================================================
            
            % y0, z0, g must be column vectors
            assert( iscolumn( y0 ) && iscolumn( z0 ) && iscolumn( g ) );

            % The length must also be identical 
            assert( length( y0 ) == length( z0 ) )
            assert( length( y0 ) == length( g  ) )

            % t_arr must be a row vector and t0i must be smaller than t_arr
            assert( isrow( t_arr ) );
            assert( t0i <= max( t_arr ) );

            % Get the number of DOF and the length of time
            n  = length( y0 );
            Nt = length( t_arr );

            % input_arr size must be just 1 smaller than Nt 
            % Since if there are N+1 steps, there should be N integration.
            [ nr, nc ] = size( input_arr );
            assert( ( n == nr ) && ( Nt == nc + 1) )
            assert( nr == length( g ) );
            
            % Define an array of y_arr, dy_arr, z_arr
            y_arr  = zeros( n, nt );
            z_arr  = zeros( n, nt );
            dy_arr = zeros( n, nt );

            % Set the initial condition 
            y_arr( :, 1 )  = y0;
            z_arr( :, 1 )  = z0;
            dy_arr( :, 1 ) = z0/obj.tau;

            for i = 1 : nt-1
                
                if t_arr( i ) <= t0i
                    y_arr( :, i+1 ) = y0;
                    z_arr( :, i+1 ) = z0;

                else
                    dt = t_arr( i+1 ) - t_arr( i );
                    [ y_new, z_new, dy, ~ ] = obj.step( y_arr( :, i ), z_arr( :, i ), g, input_arr( :, i ), dt );
                    
                    y_arr(  :, i + 1 ) = y_new;
                    z_arr(  :, i + 1 ) = z_new;
                    dy_arr( :, i + 1 ) = dy;
                end
            end

        end
    end
end

