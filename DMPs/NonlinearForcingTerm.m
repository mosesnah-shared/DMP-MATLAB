classdef NonlinearForcingTerm < handle
    
    properties

        % Canonical system and the type of Nonlinear Forcing Term
        cs
        type
        
        % Number of Basis Function of the Nonlinear Forcing Term
        N

        % The centers and widths of the N Basis Functions
        c_arr;
        h_arr;

    end

    methods
        function obj = NonlinearForcingTerm( cs, N )
            % ===========================================================================
            % Descriptions
            % ------------
            %    NonlinearForcingTerm System
            %    fs = NonlinearForcingTerm( cs, N )
            %    
            % Authors           
            % -------           
            %   Moses C. Nah    mosesnah@mit.edu
            % 
            %
            % Parameters
            % ----------
            %   (1) cs - A Canonical System
            % 
            %   (2) N  - Number of basis functions
            %
            % ===========================================================================
            
            % The canonical system of the Nonlinear Forcing term 
            obj.cs = cs; 
            
            % Type of the nonlinear forcing term is naturally defined by the Canonical System
            obj.type = obj.cs.type;
            
            % Number of Basis Functions
            % Must be positive, but if violated, MATLAB already throws an 
            % error for the next two lines 
            obj.N = N;
            obj.c_arr = zeros( 1, N );
            obj.h_arr = zeros( 1, N );

            % Setting the width and center locations of the Basis Function
            % [2023.08.15] Requires improvement for the rhythmic movement case.
            % Case for Discrete movement followed the information from:
            % Saveriano, Matteo, et al. "Dynamic movement primitives in robotics: A tutorial survey." arXiv preprint arXiv:2102.03861 (2021).

            % If Discrete Movement
            if obj.type == 0
                obj.c_arr = exp( -obj.cs.alpha_s/( obj.N-1 ) * ( 0:(obj.N-1) ) );
                obj.h_arr( 1:end-1 ) = 1.0 ./ diff( obj.c_arr ).^2;
                obj.h_arr(     end ) = obj.h_arr( end-1 );
            
            % If Rhythmic Movement
            else
                obj.c_arr = 2*pi/obj.N * ( 0.5:1:obj.N );
                obj.h_arr = obj.N * ones( 1, obj.N );
            end

        end
        
        function act = calc_ith( obj, t, i )
            % ===========================================================================
            % Descriptions
            % ------------
            %    fs = NonlinearForcingTerm.calc_ith( t, i )
            %         Calculate the activation of the i-th basis function
            %         at time t, i.e., phi_i( s( t ) )
            %    
            % Authors           
            % -------           
            %   Moses C. Nah    mosesnah@mit.edu
            % 
            %
            % Parameters
            % ----------
            %   (1) t - The time input, either scalar or array
            % 
            %   (2) i - The number of basis functions
            %
            % ===========================================================================
                       
            % The i should be within 1 to N
            assert( i >= 1 && i <= obj.N )
            
            % Getting the Center's location and width
            ci = obj.c_arr( i );
            hi = obj.h_arr( i );
            
            % The calculation of the canonical system
            s  = obj.cs.calc( t );
            
            % Calculate the activation of the i-th function
            % For discrete movement            
            if obj.type == 0
                act = Gaussian( s, ci, hi );
                
            % For rhythmic movement                
            else
                act = vonMises( s, ci, hi );
            end
        end

        function act_sum = calc_whole_at_t( obj, t )
            % ===========================================================================
            % Descriptions
            % ------------
            %    fs = NonlinearForcingTerm.calc_whole_at_t( t, i )
            %         Calculate the whole activation of the N basis
            %         functions at time t, i.e., sum_{i=1}^{N} phi_i( s( t ) )
            %    
            % Authors           
            % -------           
            %   Moses C. Nah    mosesnah@mit.edu
            % 
            %
            % Parameters
            % ----------
            %   (1) t - The time input, either scalar or array
            % ===========================================================================

            % t must be a 1D Array
            assert( length( size( t ) ) == 2 )
            [ nr, nc ] = size( t );
            assert( nr == 1 || nc == 1 ); 

            % Calculating the Activation 
            act_arr = zeros( length( t ), obj.N );

            % Iterating over the N Basis Functions
            for i = 1 : obj.N
                act_arr( :, i ) = obj.calc_ith( t, i );
            end

            % Sum over the Basis Functions
            act_sum = sum( act_arr, 2 )';

        end

        function act_weighted_sum = calc_whole_weighted_at_t( obj, t, w_arr )
            % ===========================================================================
            % Descriptions
            % ------------
            %    fs = NonlinearForcingTerm.calc_whole_at_t( t, i )
            %         Calculate the whole activation of the N basis
            %         functions at time t, i.e., sum_{i=1}^{N} phi_i( s( t ) )
            %    
            % Authors           
            % -------           
            %   Moses C. Nah    mosesnah@mit.edu
            % 
            %
            % Parameters
            % ----------
            %   (1) t - The time input, either scalar or array
            %
            %   (2) w_arr - weight array to be multiplied
            %
            % ===========================================================================

            % t must be a 1D Array
            assert( length( size( t ) ) == 2 )
            [ nr, nc ] = size( t );
            assert( nr == 1 || nc == 1 ); 

            assert( all( size( w_arr) == [ 1, obj.N ] ) );

            % Calculating the Activation 
            act_arr = zeros( length( t ), obj.N );

            % Iterating over the N Basis Functions
            for i = 1 : obj.N
                act_arr( :, i ) = w_arr( i ) * obj.calc_ith( t, i );
            end

            % Sum over the Basis Functions
            act_weighted_sum = sum( act_arr, 2 )';
        end        
    end
end