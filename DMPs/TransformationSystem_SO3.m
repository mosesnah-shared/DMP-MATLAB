classdef TransformationSystem_SO3 < handle
    
    properties
        % Parameters of the transformation systems
        alpha_z;
        beta_z;
        tau;
        
        % Initial Conditions, Orientation and Velocity
        R0;
        w0;
 
        % State values 
        R_curr;  
        w_curr;
    end
    
    methods
        function obj = TransformationSystem_SO3( alpha_z, beta_z, tau, R0, w0 )
            % ===========================================================================
            % Descriptions
            % ------------
            %    Transformation System for Orientation
            %    The orientation is represented with the R matrix, 
            %    which is an element of the SO(3) Lie Group
            %       
            %    [REF] Orientation in Cartesian space dynamic movement primitives
            %          Ude, Aleš, et al. (2014) 
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
            %   (4) R0      - Initial Orientation
            %
            %   (5) w0      - Initial Angular Velocity
            %
            % ===========================================================================
            
            assert( alpha_z > 0 && beta_z > 0 && tau > 0 );
           
            obj.alpha_z = alpha_z;
            obj.beta_z  =  beta_z;
            obj.tau     = tau;
            
            assert( obj.is_SO3( R0 ) );
            assert( isall( [ 1,3 ] == size( w0 ) ) || isall( [ 3,1 ] == size( w0 ) ) )
            
            % Initial conditions of orientation and angular velocity
            obj.R0 = R0;
            obj.w0 = w0;
            
            % Setup the current state as y0, z0 
            obj.R_curr = R0;
            obj.w_curr = w0;
            
        end
        
        % ================================================================ %
        % ===================== INTERNAL FUNCTIONS ======================= %
        % ================================================================ %
        function is_check = is_SO3( R )
            assert( isall( size( R ) == [ 3, 3 ] ) );
            is_check = ( norm( R' * R - eye( 3 ), 'fro' ) < 1e-8 );   
        end
        
        function is_check = is_so3( w_mat )
            assert( isall( size( w_mat ) == [ 3, 3 ] ) );
            is_check = ( norm( w_mat + w_mat', 'fro' ) < 1e-8 );   
        end
        
        function w_mat = R3_to_so3( w )
            assert( isall( size( w ) == [ 3, 1 ] ) && isall( size( w ) == [ 1, 3 ] ) );
            
            w_mat = [     0, -w( 3 ),  w( 2 );
                     w( 3 ),      0, -w( 1 );
                    -w( 2 ), w( 1 ),       0];     
        end
        
        function w = so3_to_R( w_mat )        
            % Should be a 3x3 matrix
            assert( is_so3( w_mat ) );
            
            w = zeros( 1, 3 );
            w( 1 ) = -w_mat( 2, 3 );
            w( 2 ) =  w_mat( 1, 3 );
            w( 3 ) = -w_mat( 1, 2 );  
        end        
        
        function R = Exp_SO3( w_mat, theta )
            % Using the Rodrigues Formula
            assert( is_so3( w_mat ) );
            
            R = eye( 3 ) + sin( theta ) * w_mat + ( 1 - cos( theta ) ) * w_mat^2;
        end
        
        function [ w_mat, theta ] = Log_so3( R )
            % First check whether the trace of R is not -1
            assert( ( trace( R ) + 1 ) >= 1e-9 );
            
            theta = acos( 0.5 * ( trace( R ) - 1 ) );
            w_mat = 1/( 2 * cos( theta ) ) * ( R - R' );
        end        
        
        % ================================================================ %
        % ================================================================ %
        % ================================================================ %        
        
        
        
        function [ R_new, w_new, dR, dw ] = step( obj, Rg, input, dt )
            % ===========================================================================            
            % Descriptions
            % ------------
            %    [ R_new, w_new, dR, d_eta ] = step( obj, Rg, input, dt )
            %       
            %    A forward integration given the following differential equation
            %       tau * dw = alpha_z * [ beta_z * log( Rg * Rc' ) - w ] + input
            %       tau * dR = [ w ] R
            %
            %    The first  equation is a dynamics on a 3D vector
            %    the second equation is simply the w represented in {S} frame
            %
            %    A simple single-step integration
            %       eta_new = eta + dt * d_eta
            %       R_new   = exp( eta_new dt/tau ) * R
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

            assert( isall( size( input ) == [1,3 ] ) || isall( size( input ) == [ 3,1 ] ) );
            
            dw = obj.alpha_z * obj.beta_z * obj.so3_to_R3( obj.Log_so3( Rg * obj.R_curr' ) ) - obj.alpha_z * obj.w + input;

            
            % Take out the y_new, z_new
            w_new = dw * dt + obj.w;
            
            w_norm = norm( w_new );
            R_new = obj.Exp_SO3( w_new/w_norm, w_norm ) * obj.R; 
            
            % Save the current R, w
            obj.R_curr = R_new;
            obj.w_curr = w_new;
            
        end
        
    end
end

