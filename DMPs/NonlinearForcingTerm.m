classdef NonlinearForcingTerm < handle
    
    properties

        % The associated Canonical system and 
        % The type of Nonlinear Forcing Term, 0: discrete, 1: rhythmic
        cs
        type  
        
        % Number of basis functions
        N

        % The centers and widths of the N basis functions
        c_arr;
        h_arr;

    end

    methods
        function obj = NonlinearForcingTerm( cs, N )
            % ===========================================================================
            % Descriptions
            % ------------
            %    Default Constructor of the Nonlinear Forcing Term
            %    fs = NonlinearForcingTerm( cs, N )
            % 
            % Parameters
            % ----------
            %   (1) cs - Canonical System
            % 
            %   (2) N  - Number of basis functions
            %
            % ===========================================================================
            
            % The Canonical System of the Nonlinear Forcing Term 
            obj.cs = cs; 
            
            % The Type of the Nonlinear Forcing Term 
            % which is naturally defined by the Canonical System
            obj.type = obj.cs.type;
            
            % Number of basis functions
            obj.N = N;
            obj.c_arr = zeros( 1, N );
            obj.h_arr = zeros( 1, N );

            % Setting the width and center locations of the Basis Function

            % For Discrete Movement, we followed from:
            % Saveriano, Matteo, et al. IJRR (2023)
            if obj.cs.type == 0
                obj.c_arr = exp( -obj.cs.alpha_s/( obj.N-1 ) * ( 0:(obj.N-1) ) );
                obj.h_arr( 1:end-1 ) = 1.0 ./ diff( obj.c_arr ).^2;
                obj.h_arr(     end ) = obj.h_arr( end-1 );
            
            % For Rhythmic Movement, the ci array followed from:
            % Saveriano, Matteo, et al. IJRR (2023)            
            else
                obj.c_arr = 2*pi/obj.N * ( 0.5:1:obj.N );
                obj.h_arr = obj.N * ones( 1, obj.N );
            end

        end
        
        function act = calc_ith( obj, t_arr, i )
            % ===========================================================================
            % Descriptions
            % ------------
            %    fs = NonlinearForcingTerm.calc_ith( t_arr, i )
            %         Calculate the activation of the i-th basis function
            %         at time t_arr (not s_arr of Canonical System)
            %         Meaning, phi_i( s( t_arr ) )
            % 
            % Parameters
            % ----------
            %   (1) t_arr - The time array row vector.
            % 
            %   (2) i     - The i-th basis functions, an integer.
            %
            % Returns
            % -------
            %   (1) act   - The activation array of phi_i
            %
            % ===========================================================================

            % t must be a row vector
            assert( isrow( t_arr ) )     

            % The i must be scalar value within 1 to N
            assert( isscalar( i ) && i >= 1 && i <= obj.N )
            
            % Getting the center location and width of the basis function
            ci = obj.c_arr( i );
            hi = obj.h_arr( i );
            
            % The calculation of the Canonical system
            s_arr  = obj.cs.calc( t_arr );
            
            % Calculate the activation of the i-th function
            % For discrete movement            
            if obj.type == 0
                act = Gaussian( s_arr, ci, hi );
                
            % For rhythmic movement                
            elseif obj.type == 1
                act = vonMises( s_arr, ci, hi );

            else
                error( 'Wrong type for the Nonlinear Forcing Term' )
            end
        end

        function act_arr = calc_multiple_ith( obj, t_arr, i_arr )
            % ===========================================================================
            % Descriptions
            % ------------
            %    fs = NonlinearForcingTerm.calc_multiple_ith( t_arr, i_arr )
            %         Calculate the activation of multiple i-th 
            %         basis functions at time t.
            %         Internally, it called calc_ith method multiple times.
            %    
            % 
            % Parameters
            % ----------
            %   (1) t_arr - The time input array
            % 
            %   (2) i_arr - The array of basis functions
            %
            % Returns
            % -------
            %   (1) act_arr - The activation array of multiple phi_i
            %                 Size N x Nt, where:
            %                 (1)  N: The number of basis functions
            %                 (2) Nt: The length of t_arr
            %
            % ===========================================================================
                       
            % Both i_arr and t_arr must be row vector
            assert( isrow( i_arr ) && isrow( t_arr ), 'Both input must be row arrays' )         

            % Sorting the i_arr in case if not sorted
            i_arr = sort( i_arr );
            
            % Initialize the activation array 
            %    row: # of i-array
            % column: # of time-array
            act_arr = zeros( length( i_arr ), length( t_arr ) );
            
            % Iterate along the basis functions.
            for i = 1 : length( i_arr )
                act_arr( i, : ) = obj.calc_ith( t_arr , i_arr( i ) );
            end

        end        


        function act_sum = calc_whole_at_t( obj, t_arr )
            % ===========================================================================
            % Descriptions
            % ------------
            %    fs = NonlinearForcingTerm.calc_whole_at_t( t_arr )
            %         Calculate the whole activation of the N basis
            %         functions at time t, i.e., sum_{i=1}^{N} phi_i( s( t_arr ) )
            %         This is the denominator of the Nonlinear Forcing Term
            %
            % Parameters
            % ----------
            %   (1) t_arr - The time array input
            %
            % Returns
            % -------          
            %   (1) act_sum - The array of sum of activation
            % ===========================================================================

            % Sum over the activation of N basis functions
            % Sum along the row, which is along basis functions
            act_sum = sum( obj.calc_multiple_ith( t_arr, 1:obj.N ), 1 );

        end

        function act_weighted = calc_whole_weighted_at_t( obj, t_arr, w_arr )
            % ===========================================================================
            % Descriptions
            % ------------
            %    fs = NonlinearForcingTerm.calc_whole_weighted_at_t( t_arr, w_arr )
            %         Calculate the whole activation of the N basis
            %         with the weighting array w_arr
            %         This is the numerator of the Nonlinear Forcing Term            
            % 
            % Parameters
            % ----------
            %   (1) t_arr - The time row array input
            %
            %   (2) w_arr - The 2D Weight array to be multiplied
            %               n x N, where 
            %               (1) n: Number of dimension of the Nonlinear Forcing Term
            %               (2) N: Number of basis functions.
            %
            % Returns
            % -------
            %   (1) act_weighted - sum_{i=1}^{N}w_i phi_i( s( t ) )/sum_{i=1}^{N}phi_i( s( t ) )
            %       
            % ===========================================================================

            % t must be a row vector
            assert( isrow( t_arr ) )
            
            % The # of column of the Weight array must be identical 
            % to the number of basis functions
            [ n, nc ] = size( w_arr );
            assert( nc == obj.N );
            
            % Initialize the Nonlinear Forcing Term 
            act_weighted = zeros( n, length( t_arr ) );

            for i = 1 : n
                act_weighted( i, : ) = sum( w_arr( i, : )' .* obj.calc_multiple_ith( t_arr, 1:obj.N ) ./obj.calc_whole_at_t( t_arr ), 1 );
            end
        end   

        function force_arr = calc_forcing_term( obj, t_arr, w_arr, t0i, scl )
            % ===========================================================================
            % Descriptions
            % ------------
            %    fs = NonlinearForcingTerm.calc_forcing_term( t_arr, w_arr, t0i, scl )
            %         Calculate the nonlinear forcing term with a time offset of t0i
            %         This function is for Imitation Learning
            %
            % Parameters
            % ----------
            %   (1) t_arr - The time row array input
            %
            %   (2) w_arr - The 2D Weight array to be multiplied
            %               n x N, where 
            %               (1) n: Number of dimension of the Nonlinear Forcing Term
            %               (2) N: Number of basis functions.
            %
            %   (3) t0i   - The time offset for the nonlinear forcing term
            %
            %   (4) scl   - The 2D scaling square matrix
            %               Must be the size of n x n
            %
            % Returns
            % -------
            %   (1) force_arr - The n x Nt, where 
            %                   Nt is the length of the time array.
            %     
            % ===========================================================================
    
            % t must be a row vector, with length Nt
            assert( isrow( t_arr ) )        
            Nt = length( t_arr );
    
            % Number of column of weight array must be equal to 
            % the number of basis functions N
            [ n, nc ] = size( w_arr );
            assert( nc == obj.N );

            % Check the scaling matrix is nxn matrix
            assert( all( [ n, n ] == size( scl ) ) );
                
            % The time array's maximum value should be bigger than t0i
            assert( t0i <= max( t_arr ) && t0i >= 0 );
    
            % Initialize the forcing term array
            force_arr = zeros( n, Nt );
    
            % Get the index for the time array
            % That is between t0i and t0i + D
            %idx_arr = ( t_arr >= t0i & t_arr <= t0i + obj.cs.tau );
            idx_arr = ( t_arr >= t0i  );

            % Shifting the element one side to the right 
            % Adding this is quite important, since 
            % The rollout requires one-step further from the integration
            % [Notes] [Moses C. Nah] [2023.11.16]
            idx_arr = circshift( idx_arr, -1 );
    
            % Calculate the forcing term 
            t_tmp = t_arr( idx_arr ) - t0i;
            force_arr( :, idx_arr ) = scl * ( obj.cs.calc( t_tmp ) .* obj.calc_whole_weighted_at_t( t_tmp, w_arr ) );

        end

    end
end