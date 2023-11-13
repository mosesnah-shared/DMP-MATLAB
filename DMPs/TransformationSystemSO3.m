classdef TransformationSystemSO3 < handle
    
    properties
        % Parameters of the transformation systems
        alpha_z;
        beta_z;
        tau;
        
        % The Canonical System, to get the tau value
        % Note that the Nonlinear Forcing Term is not associated with the
        % Transformation System
        cs; 

    end
    
    methods
        function obj = TransformationSystemSO3( alpha_z, beta_z, cs )
            % ===========================================================================
            % Descriptions
            %    Constructor of Transformation System
            %    This transformation system is for Orientation with
            %    SO(3) Representation
            %
            %    trans_sys = TransformationSystem( alpha_z, beta_z, cs )
            %       
            %    It is a second-order linear system with an input
            %    as a Nonlinear Forcing Term
            %       tau * dR( t ) = [ w(t) ]R( t )
            %       tau * dw( t ) = alpha_z * [ beta_z * Log( R ) - w( t ) ] + input
            %
            % Parameters
            % ----------
            %   (1) alpha_z - A Positive Gain
            % 
            %   (2) beta_z  - A Positive Gain
            %
            %   (3) cs      - Canonical System
            %
            % ===========================================================================
            
            assert( alpha_z > 0 && beta_z > 0 );
           
            obj.alpha_z = alpha_z;
            obj.beta_z  =  beta_z;
            obj.cs      =  cs;
            obj.tau     =  obj.cs.tau;
            
        end
       
        
        function [ R_new, w_new, dR, dw ] = step( obj, R_old, w_old, input, dt )
            % ===========================================================================            
            % Descriptions
            % ------------
            %    [ R_new, w_new, dR, dw ] = step( obj, R_old, w_old, input, dt )
            %       
            %    A forward integration given the following differential equation
            %       tau * dR = [ w_old ] R
            %       tau * dw = alpha_z * ( beta_z ( Log( Rold ) - w ) + input 
            %
            % Parameters
            % ----------
            %   (1) R_old - The SO3 matrix
            % 
            %   (2) w_old - The associated angular velocity
            %
            %   (3) input - The input forcing term 
            %
            %   (4) dt    - time step
            %
            % Returns
            % -------
            %   (1) R_new - The integrated SO3 matrix
            %
            %   (2) w_new - The integrated angular velocity 
            %
            %   (3) dR    - time derivative of R
            %
            %   (4) dw    - time derivative of angular velocity
            %
            % ===========================================================================

            % Check whether the input and w_old are column vectors
            assert( iscolumn( input ) && iscolumn( w_old ) );
            assert(   length( input ) ==   length( w_old ) );

            % The SO3 matrix
            assert( is_SO3( R_old ) );

            % Time derivative of R and dw
            dw = 1/obj.tau * obj.alpha_z * ( obj.beta_z * so3_to_R3( LogSO3( R_old ) ) - w_old ) + input;
                
            % The forward integration
            R_new = dR * dt + R_old; % Use
            w_new = dw * dt + w_old;

        end

        function f_des = get_desired( obj, R_arr, w0_arr, dw0_arr )
            % Descriptions
            % ------------
            %    Calculate the Nonlinear Forcing Terms that should be
            %    learned for Imitation Learning
            %
            % Parameters
            % ----------
            %   (1) R_arr   - An array of SO3 matrices
            %                 3x3xNp of data
            % 
            %   (2) w0_arr  - Angular velocity, 3xNp
            %
            %   (3) dw0_arr - Angular acceleration, 3xNp
            %
            % Returns
            % -------
            %   (1) f_des = tau * dw0_arr + alpha_z * w0_arr + alpha_z beta_z Log (R_arr)
            %
            % ===========================================================================

            % Check the length of the w0_arr, dw0_arr and R_arr
            [ ~, ~, Np  ] = size(  R_arr  );
            [ ~,    Np1 ] = size(  w0_arr );
            [ ~,    Np2 ] = size( dw0_arr );

            assert( Np == Np1 && Np == Np2 )

            % Get the Log of R_arr 
            logR = zeros( 3, Np );
            for i = 1 : Np
                R = R_arr( :, :, i );
                logR( :, i ) = so3_to_R3( LogSO3( R ) );
            end
            f_des = obj.tau * dw0_arr + obj.alpha_z * w0_arr + obj.alpha_z * obj.beta_z * logR;
            
        end
    end
end

