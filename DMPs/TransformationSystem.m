classdef TransformationSystem < handle
    
    properties
        alpha_z;
        beta_z;
        
        tau;
        y0;
        z0;
 
        % State values 
        y_curr;
        z_curr;
    end
    
    methods
        function obj = TransformationSystem( alpha_z, beta_z, tau, y0, z0 )
            % Setting up the transformation system
            obj.alpha_z = alpha_z;
            obj.beta_z  =  beta_z;
            
            obj.y0 = y0;
            obj.z0 = z0;
            
            obj.tau = tau;
            
            % Setup the current state as y0, z0 
            obj.y_curr = y0;
            obj.z_curr = z0;
            
        end
        
        function [ y_new, z_new, dy, dz ] = step( obj, g, input, dt )
            % Forward integrate with the goal value, input value 
            dy = obj.z_curr / obj.tau;
            dz = ( obj.alpha_z * obj.beta_z * ( g - obj.y_curr ) - obj.alpha_z * obj.z_curr + input ) / obj.tau;
            
            % Take out the y_new, z_new
            y_new = dy * dt + obj.y_curr;
            z_new = dz * dt + obj.z_curr;
            
            obj.y_curr = y_new;
            obj.z_curr = z_new;
            
            
        end
    end
end

