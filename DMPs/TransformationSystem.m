classdef TransformationSystem < handle
    
    properties
        
        % Parameters of the transformation systems
        alpha_z;
        beta_z;
        tau;
        
        % Initial Conditions
        y0;
        z0;
 
        % State values 
        y_curr;
        z_curr;
    end
    
    methods
        function obj = TransformationSystem( alpha_z, beta_z, tau, y0, z0 )
            % ===========================================================================
            % Descriptions
            % ------------
            %    Transformation System
            %    trans_sys = TransformationSystem( alpha_z, beta_z, tau, y0, z0 )
            %       
            %    It is a second-order linear system with an input
            %       tau * dy = z
            %       tau * dz = alpha_z * [ beta_z * (g - y) - z ] + input
            %    
            %    Basic Implementation of Transformation System
            %
            % Authors           
            % -------           
            %   Moses C. Nah    mosesnah@mit.edu
            % 
            %
            % Parameters
            % ----------
            %   (1) alpha_z - A Positive Gain
            % 
            %   (2) beta_z  - A Positive Gain
            %
            %   (3) tau     - Time constant, derived from the Canonical System
            %
            %   (4) y0      - Initial position
            %
            %   (5) z0      - Initial (normalized) velocity
            %
            % ===========================================================================
            
            assert( alpha_z > 0 && beta_z > 0 && tau > 0 );
           
            obj.alpha_z = alpha_z;
            obj.beta_z  =  beta_z;
            obj.tau     = tau;
            
            % Initial conditions
            obj.y0 = y0;
            obj.z0 = z0;
            
            % Setup the current state as y0, z0 
            obj.y_curr = y0;
            obj.z_curr = z0;
            
        end
        
        function [ y_new, z_new, dy, dz ] = step( obj, g, input, dt )
            % ===========================================================================            
            % Descriptions
            % ------------
            %    [ y_new, z_new, dy, dz ] = step( obj, g, input, dt )
            %       
            %    A forward integration given the following differential equation
            %       tau * dy = z
            %       tau * dz = alpha_z * [ beta_z * (g - y) - z ] + input
            %    
            %    A simple single-step integration
            %       y_new = y + dt * dy   
            %       z_new = z + dt * dz   
            % 
            % Authors           
            % -------           
            %   Moses C. Nah    mosesnah@mit.edu
            % 
            %
            % Parameters
            % ----------
            %   (1) g     - Goal position, could be modified
            % 
            %   (2) input - force, could include coupling term
            %
            %   (3) dt    - time step
            %
            % ===========================================================================

            dy = obj.z_curr/obj.tau;
            dz = ( obj.alpha_z * obj.beta_z * ( g - obj.y_curr ) - obj.alpha_z * obj.z_curr + input )/obj.tau;
            
            % Take out the y_new, z_new
            y_new = dy * dt + obj.y_curr;
            z_new = dz * dt + obj.z_curr;
            
            % Save the current y,z
            obj.y_curr = y_new;
            obj.z_curr = z_new;
            
        end
        
        function reset( obj )
            % Reset the current state as y0, z0 
            obj.y_curr = obj.y0;
            obj.z_curr = obj.z0;
        end
    
        
        function f_des = get_desired( obj, y_des, dy_des, ddy_des, g )
            % ===========================================================================            
            % Descriptions
            % ------------
            %    f_des = get_desired( obj, y_des, dy_des, ddy_des, g )
            %       
            %    Calculate an Array of f_des, from the y_desired Traj.
            % 
            % Authors           
            % -------           
            %   Moses C. Nah    mosesnah@mit.edu
            % 
            % Parameters
            % ----------
            %   (1)   y_des -     Position of Desired Trajectory, 1xP array
            % 
            %   (2)  dy_des -     Velocity of Desired Trajectory, 1xP array
            %
            %   (3) ddy_des - Acceleration of Desired Trajectory, 1xP array
            %
            %   (4)       g - goal location, a constant
            % ===========================================================================
                
            % Check whether y_des, dy_des, ddy_des are all same size 
            assert( all( size( y_des ) == size(  dy_des ) ) );
            assert( all( size( y_des ) == size( ddy_des ) ) );

            % Check that y_des are 1D array
            [nr, nc] = size( y_des );
            assert( nr == 1 || nc == 1);

            % Get the length of desired trajectory, i.e., 
            % The number of sample points P
            P = length( y_des );

            % Calculate the f_desired array
            f_des = obj.tau^2 * ddy_des + obj.alpha_z * obj.tau * dy_des + obj.alpha_z * obj.beta_z * ( y_des - g );

        end
    end
end

