classdef TransformationSystemQuat < handle
    
    properties
        % Parameters of the transformation systems
        alpha_z;
        beta_z;
        tau;
        
        % Initial Conditions, Orientation and Velocity
        quat0;
        w0;
 
        % State values 
        quat_curr;  
        w_curr;
    end
    
    methods
        function obj = TransformationSystemQuat( alpha_z, beta_z, tau, quat0, w0 )
            % ===========================================================================
            % Descriptions
            % ------------
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
            %   (4) quat0   - Initial Orientation of the Quaternion
            %
            %   (5) w0      - Initial Angular Velocity
            %
            % ===========================================================================
            
            assert( alpha_z > 0 && beta_z > 0 && tau > 0 );
           
            obj.alpha_z = alpha_z;
            obj.beta_z  =  beta_z;
            obj.tau     = tau;
            
            assert( obj.is_unit_quat( quat0 ) );
            assert( all( [ 1,3 ] == size( w0 ) ) || all( [ 3,1 ] == size( w0 ) ) )
            
            if all( [ 1,3 ] == size( w0 )  )
               w0 = w0'; 
            end
            
            if all( [ 1,4 ] == size( quat0 )  )
               quat0 = quat0'; 
            end
                        
            % Initial conditions of orientation and angular velocity
            obj.quat0 = quat0;
            obj.w0    = w0;
            
            % Setup the current state as y0, z0 
            obj.quat_curr = quat0;
            obj.w_curr    =    w0;
            
        end
        
        % ================================================================ %
        % ===================== INTERNAL FUNCTIONS ======================= %
        % ================================================================ %
        function is_check = is_unit_quat( obj, quat )
            assert( all( size( quat ) == [ 1, 4 ] ) || all( size( quat ) == [ 4, 1 ] ) );
            assert( abs( dot( quat, quat ) - 1 ) <= 1e-9 );
            is_check = true;
        end
        
        
        function [ quat_new, w_new, dquat, dw ] = step( obj, quatg, input, dt )
            % ===========================================================================            
            % Descriptions
            % ------------
            %    [ quat_new, w_new, dquat, dw ] = step( obj, quatg, input, dt )
            %       
            %    A forward integration given the following differential equation
            %       tau * dq = 1/2 w * q
            %       tau * dw = alpha_z * ( beta_z ( vec( quatg * q' ) - w )
            %
            %    The first  equation is the Quaternion Dynamics
            %    w in this is simply for w = (wx, wy, wz):
            %    w = wxi + wyj + wzk
            %    the second equation is simply the w.
            %    vec is taking off the imaginary part
            %
            %    [REF1] Saveriano, Matteo, et al. "Dynamic movement primitives in robotics: A tutorial survey." 
            %           arXiv preprint arXiv:2102.03861 (2021).
            %    [REF2] Saveriano, Matteo, Felix Franzel, and Dongheui Lee. 
            %           "Merging position and orientation motion primitives." (2019) Eq.3b
            % 
            % Authors           
            % -------          
            %   Moses C. Nah    mosesnah@mit.edu
            % 
            %
            % Parameters
            % ----------
            %   (1) Rg    - Goal orientation of the robot
            % 
            %   (2) input - force, could include coupling term
            %
            %   (3) dt    - time step
            %
            % ===========================================================================

            assert( all( size( input ) == [ 1,3 ] ) || all( size( input ) == [ 3,1 ] ) );
            assert( obj.is_unit_quat( quatg ) );
                
            dquat      = 1/2 * quat_mul( [ 0; obj.w_curr ], obj.quat_curr ); 
            quat_delta = quat_log( quat_mul( quatg, quat_conj( obj.quat_curr ) ) );
            
            dw = ( obj.alpha_z * obj.beta_z * quat_delta - obj.beta_z * obj.w_curr + input )/obj.tau;
            
            w_new  = dw * dt + obj.w_curr;
            
            % Forward integration of Quaternion
            quat_new = quat_mul( quat_exp( 1/2 * dt/obj.tau * obj.w_curr ), obj.quat_curr );
            
            % Update the current one
            obj.quat_curr = quat_new;
            obj.w_curr    = w_new;
           

        end
        
    end
end

