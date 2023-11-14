classdef TransformationSystemQuat < handle
    
    properties

        % Parameters of the Transformation System
        % The gains are all identical over the 3-DOFs
        alpha_z;
        beta_z;
        tau;

        % The Canonical System, to get the tau value
        % Note that the Nonlinear Forcing Term is not associated with the
        % Transformation System
        cs; 
 
    end
    
    methods
        function obj = TransformationSystemQuat( alpha_z, beta_z, cs )
            % ===========================================================================
            % Descriptions
            % ------------
            %    Constructor of Transformation System for Orientation, Unit Quaternion description 
            %    The method is derived by Koutras and Doulgeri (2020)
            %    "A correct formulation for the orientation dynamic movement primitives for robot control in the cartesian space." 
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
            
            % Must be positive values
            assert( alpha_z > 0 && beta_z > 0 );
           
            obj.alpha_z = alpha_z;
            obj.beta_z  =  beta_z;
            obj.cs      = cs;
            obj.tau     = cs.tau;
            
        end
                
        function [ eq_new, z_new, deq, dz ] = step( obj, eq_old, z_old, input, dt )
            % ===========================================================================            
            % Descriptions
            % ------------
            %    [ eq_new, z_new, deq_new, dz ] = step( obj, eq_old, z_old, input, dt )
            %       
            %    A forward integration given the following differential equation
            %       tau * dz  = eq_old
            %       tau * deq = -alpha_z * ( beta_z eq_old + z ) + input
            %
            % Parameters
            % ----------
            %   (1) eq_old - error quat.
            %             
            %   (2) z_old  - the generalized time derivative of error quat.
            % 
            %   (3) input  - Nonlinear Forcing Term input
            %
            %   (4) dt     - Time step
            %
            % Returns
            % -------
            %   (1) eq_new - Integrated new position column array
            % 
            %   (2)  z_new - Integrated new generalized velocity array
            %
            %   (3)    deq - The velocity column array, z/tau
            %
            %   (4)     dz - The acceleration column array
            %
            % ===========================================================================

            assert( iscolumn( eq_old ) && iscolumn( z_old ) && iscolumn( input ) )
            assert( ( length( eq_old ) == 3 ) && ( length( z_old ) == 3 ) && ( length( input ) == 3  ) )
                
            deq = 1/obj.tau * z_old;
            dz  = 1/obj.tau * ( -obj.alpha_z * ( obj.beta_z * eq_old + z_old ) + input );

            z_new  =  z_old + dz  * dt;
            eq_new = eq_old + deq * dt;

        end


        function [ eq_arr, z_arr, deq_arr ] = rollout( obj, eq0, z0, input_arr, t0i, t_arr )
            % ===========================================================================            
            % Descriptions
            % ------------
            %    [eq_arr, z_arr, deq_arr ] = rollout( obj, eq0, deq0, input_arr, t0i, t_arr )
            %       
            %    A forward integration given the following differential equation
            %       tau * dz  = eq_old
            %       tau * deq = -alpha_z * ( beta_z eq_old + z ) + input
            %
            % Parameters
            % ----------
            %   (1) eq0  - The initial condition of error quaternion
            %             
            %   (2) z0   - The initial condition of z
            %
            %   (3) input_arr - Nonlinear Forcing Term input array
            %
            %   (4) t0i - Initial time of the rollout
            %             Must be identical with input_arr
            % 
            %   (5) t_arr - time array of the integration
            %
            % Returns
            % -------
            %   (1) eq_arr  - Error quaternion array
            % 
            %   (2) z_arr   - deq_arr / tau
            %
            %   (3) deq_arr - Time-derivative of Error quaternion array
            %
            % ===========================================================================            
            % eq0, z0 must be column vectors
            assert( iscolumn( eq0 ) && iscolumn( z0 ) );

            % The length must also be identical 
            assert( length( eq0 ) == length( z0 ) )

            % t_arr must be a row vector and 
            % t0i must be smaller than the maximum value of t_arr 
            assert( isrow( t_arr ) );
            assert( isscalar( t0i ) && t0i <= max( t_arr ) );

            % The length of time
            Nt = length( t_arr );

            % input_arr size must be just 1 smaller than Nt 
            % Since if there are Nt steps, there should be Nt-1 integrations
            [ nr, nc ] = size( input_arr );
            assert( ( 3 == nr ) && ( Nt == nc + 1 ) )
            
            % Initialize eq_arr, z_arr, deq_arr for the integration
            eq_arr  = zeros( 3, Nt );
            z_arr   = zeros( 3, Nt );
            deq_arr = zeros( 3, Nt );

            % Set the initial condition 
            eq_arr( :, 1 )  = eq0;
            z_arr( :, 1 )   = z0;
            deq_arr( :, 1 ) = z0/obj.tau;

            for i = 1 : Nt-1
                
                % If smaller than t0i, no integration
                if t_arr( i+1 ) <= t0i
                    eq_arr( :, i+1 ) = eq0;
                    z_arr( :, i+1 )  = z0;

                else
                    dt = t_arr( i+1 ) - t_arr( i );
                    [ y_new, z_new, dy, ~ ] = obj.step( eq_arr( :, i ), z_arr( :, i ), input_arr( :, i ), dt );
                    
                    eq_arr(  :, i + 1 ) = y_new;
                    z_arr(  :, i + 1 )  = z_new;
                    deq_arr( :, i + 1 ) = dy;
                end
            end            

        end        

        function f_des = get_desired( obj, eq_arr, deq_arr, ddeq_arr )
            % ===========================================================================            
            % Descriptions
            % ------------
            %    f_des = get_desired( obj, eq_arr, deq_arr, ddeq_arr )
            %       
            %    Calculating the f_des for Imitation Learning
            %
            % Parameters
            % ----------
            %   (1) eq_arr   - The error quaternion array, 3xNq
            %             
            %   (2) deq_arr  - The time-derivative  of error quaternion array, 3xNq
            %             
            %   (3) ddeq_arr - The time-dderivative of error quaternion array, 3xNq
            %
            % Returns
            % -------
            %   (1) f_des - The desired force for the quaternion transformation system
            %
            % ===========================================================================            
                        

            f_des = obj.tau^2 * ddeq_arr + obj.alpha_z * obj.tau * deq_arr + obj.alpha_z * obj.beta_z * eq_arr;

        end

        
    end
end

