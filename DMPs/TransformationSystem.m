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
        
        
    end
end

