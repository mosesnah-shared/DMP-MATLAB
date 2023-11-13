classdef TransformationSystem < handle
    
    properties
        
        % Parameters of the Transformation System
        % The gains are all identical over the n-DOF
        alpha_z;
        beta_z;
        tau;

        % The Canonical System, to get the tau value
        % Note that the Nonlinear Forcing Term is not associated with the
        % Transformation System
        cs; 
 
    end
    
    methods
        function obj = TransformationSystem( alpha_z, beta_z, cs )
            % ===========================================================================
            % Descriptions
            % ------------
            %    Constructor of Transformation System
            %    This transformation system is for joint-space and task-space, position.    
            %
            %    trans_sys = TransformationSystem( alpha_z, beta_z, cs )
            %       
            %    It is a second-order linear system with an input
            %    as a Nonlinear Forcing Term
            %       tau * dy( t ) = z( t )
            %       tau * dz( t ) = alpha_z * [ beta_z * ( g - y( t ) ) - z( t ) ] + input
            % 
            % Parameters
            % ----------
            %   (1) alpha_z - A Positive Gain
            % 
            %   (2) beta_z  - A Positive Gain
            %
            %   (3) cs      - Canonical System, 
            %                 which provides the tau value of the system
            %
            % ===========================================================================
            
            % Must be positive values
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
            %    A forward integration of the following differential equation
            %       tau * dy = z_old
            %       tau * dz = alpha_z * [ beta_z * ( g - y_old ) - z_old ] + input
            %    
            %    A single-step integration with time-step dt
            %       y_new = y_old + dt * dy   
            %       z_new = z_old + dt * dz   
            %
            % Parameters
            % ----------
            %   (1) y_old - The position column array
            %             
            %   (2) z_old - The generalized velocity column array
            %             
            %   (3) g     - The goal position column array
            %               Note that the goal position can be 
            %               a function of time
            % 
            %   (4) input - Nonlinear Forcing Term input
            %
            %   (5) dt    - Time step
            %
            % Returns
            % -------
            %   (1) y_new - Integrated new position column array
            % 
            %   (2) z_new - Integrated new generalized velocity array
            %
            %   (3) dy    - The velocity column array, z/tau
            %
            %   (4) dz    - The acceleration column array
            %
            % ===========================================================================
            
            % For the Transformation Systems,
            % The arguments must be all column arrays
            assert( iscolumn( y_old ) && iscolumn( z_old ) && ...
                    iscolumn(     g ) && iscolumn( input ) )
            
            % The length must all be identical 
            assert( length( y_old ) == length( z_old )  )
            assert( length( y_old ) == length( g     )  )
            assert( length( y_old ) == length( input )  )
            
            % Time-step must be strictly positive 
            assert( isscalar( dt ) && dt > 0 )
         
            % Time-derivative of y, z 
            dy = z_old/obj.tau;
            dz = 1/obj.tau * ( obj.alpha_z * obj.beta_z * ( g - y_old ) - obj.alpha_z * z_old + input );
            
            % Forward Integration
            y_new = dy * dt + y_old;
            z_new = dz * dt + z_old;
            
        end
        
        function [ y_arr, z_arr, dy_arr ] = rollout( obj, y0, z0, g, input_arr, t0i, t_arr )
            % ===========================================================================            
            % Descriptions
            % ------------
            %    [ y_arr, z_arr, dy_arr ] = rollout( obj, y0, z0, g, input_arr, t0i, t_arr )
            %       
            %    A full rollout of the trajectory
            %       tau * dy = z_old
            %       tau * dz = alpha_z * [ beta_z * ( g - y_old ) - z_old ] + input
            %    
            %    The forward integration
            %       y_new = y_old + dt * dy   
            %       z_new = z_old + dt * dz   
            %
            % Parameters
            % ----------
            %   (1) y0 - The initial condition of the position, column array
            %             
            %   (2) z0 - The initial condition of the generalized velocity column array
            %             
            %   (3) g  - The goal position column array
            %            Constant vector.
            % 
            %   (4) input_arr - Nonlinear Forcing Term input array
            %
            %   (5) t0i - Initial time of the rollout
            %             Must be identical with input_arr
            % 
            %   (6) t_arr - time array of the integration
            %
            % Returns
            % -------
            %   (1) y_arr  - The position, concatenated along column array
            % 
            %   (2) z_arr  - The generalized velocity, concatenated along column array
            %
            %   (3) dy_arr - The velocity, concatenated along column array
            %
            % ===========================================================================            
            
            % y0, z0, g must be column vectors
            assert( iscolumn( y0 ) && iscolumn( z0 ) && iscolumn( g ) );

            % The length must also be identical 
            assert( length( y0 ) == length( z0 ) )
            assert( length( y0 ) == length( g  ) )

            % t_arr must be a row vector and 
            % t0i must be smaller than the maximum value of t_arr 
            assert( isrow( t_arr ) );
            assert( isscalar( t0i ) && t0i <= max( t_arr ) );

            % Get the number of DOF and the length of time
            n  = length( y0 );
            Nt = length( t_arr );

            % input_arr size must be just 1 smaller than Nt 
            % Since if there are Nt steps, there should be Nt-1 integrations
            [ nr, nc ] = size( input_arr );
            assert( ( n == nr ) && ( Nt == nc + 1 ) )
            
            % Initialize y_arr, dy_arr, z_arr for the integration
            y_arr  = zeros( n, Nt );
            z_arr  = zeros( n, Nt );
            dy_arr = zeros( n, Nt );

            % Set the initial condition 
            y_arr( :, 1 )  = y0;
            z_arr( :, 1 )  = z0;
            dy_arr( :, 1 ) = z0/obj.tau;

            for i = 1 : Nt-1
                
                % If smaller than t0i, no integration
                if t_arr( i+1 ) <= t0i
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

        function f_des = get_desired( obj, y_des, dy_des, ddy_des, g )
            % ===========================================================================            
            % Descriptions
            % ------------
            %    f_des = get_desired( obj, y_des, dy_des, ddy_des, g )
            %       
            %    Calculate an array of f_des for Imitation Learning
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
            %
            % Returns
            % -------
            %   (1) f_des   - tau^2 ddy_des + alpha_z dy_des + alpha_z beta_z ( y_des - g )
            %               - n x P array
            %
            % ===========================================================================
                
            % Check whether the 2D arrays of y_des, dy_des, ddy_des are all same size 
            assert( all( size( y_des ) == size(  dy_des ) ) );
            assert( all( size( y_des ) == size( ddy_des ) ) );
            assert( iscolumn( g ) )

            % Check whether the size of goal is consistent 
            [ n, ~ ] = size( y_des );
            assert( length( g ) == n );

            % Calculate the f_des array
            f_des = obj.tau^2 * ddy_des + obj.alpha_z * obj.tau * dy_des + obj.alpha_z * obj.beta_z * ( y_des - g );

        end        

    end
end

