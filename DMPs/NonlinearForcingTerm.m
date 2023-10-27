classdef NonlinearForcingTerm < handle
    
    properties

        % Canonical system and the type of Nonlinear Forcing Term
        cs
        type  % 0: discrete, 1: rhythmic
        
        % Number of basis functions of the Nonlinear Forcing Term
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
            %    NonlinearForcingTerm System
            %    fs = NonlinearForcingTerm( cs, N )
            %    
            % Authors           
            % -------           
            %   Moses C. Nah    mosesnah@mit.edu
            % 
            % Parameters
            % ----------
            %   (1) cs - Canonical System
            % 
            %   (2) N  - Number of basis functions
            %
            % ===========================================================================
            
            % The canonical system of the Nonlinear Forcing term 
            obj.cs = cs; 
            
            % Type of the nonlinear forcing term is naturally defined by the Canonical System
            obj.type = obj.cs.type;
            
            % Number of basis functions
            obj.N = N;
            obj.c_arr = zeros( 1, N );
            obj.h_arr = zeros( 1, N );

            % Setting the width and center locations of the Basis Function
            % [2023.08.15] Requires improvement for the rhythmic movement case.
            % Case for Discrete movement followed the information from:
            % Saveriano, Matteo, et al. "Dynamic movement primitives in robotics: A tutorial survey." arXiv preprint arXiv:2102.03861 (2021).

            % If discrete Movement
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
        
        function act = calc_ith( obj, t_arr, i )
            % ===========================================================================
            % Descriptions
            % ------------
            %    fs = NonlinearForcingTerm.calc_ith( t_arr, i )
            %         Calculate the activation of the i-th basis function
            %         at time t_arr, i.e., phi_i( s( t_arr ) )
            %    
            % Authors           
            % -------           
            %   Moses C. Nah    mosesnah@mit.edu
            % 
            % Parameters
            % ----------
            %   (1) t_arr - The time input row vector.
            % 
            %   (2) i - The number of basis functions, an integer.
            %
            % Returns
            % -------
            %   (1) act - The activation array of phi_i
            %
            % ===========================================================================
                       
            % The i should be within 1 to N
            assert( i >= 1 && i <= obj.N )

            % t must be a row vector
            assert( isrow( t_arr ) )     
            
            % Getting the Center's location and width
            ci = obj.c_arr( i );
            hi = obj.h_arr( i );
            
            % The calculation of the canonical system
            s_arr  = obj.cs.calc( t_arr );
            
            % Calculate the activation of the i-th function
            % For discrete movement            
            if obj.type == 0
                act = Gaussian( s_arr, ci, hi );
                
            % For rhythmic movement                
            else
                act = vonMises( s_arr, ci, hi );
            end
        end

        function act_arr = calc_ith_arr( obj, t_arr, i_arr )
            % ===========================================================================
            % Descriptions
            % ------------
            %    fs = NonlinearForcingTerm.calc_ith_arr( t_arr, i_arr )
            %         Calculate the activation of the i_arrs basis function
            %         basis functions at time t, i.e., phi_i( s( t ) )
            %    
            % Authors           
            % -------           
            %   Moses C. Nah    mosesnah@mit.edu
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
            %                 Size N x Nt, where Nt is the length of t_arr
            %
            % ===========================================================================
                       
            % Both i_arr and t_arr must be row vector
            assert( isrow( i_arr ) && isrow( t_arr ) )         

            % Sorting the i_arr 
            i_arr = sort( i_arr );
            
            % Initialize the activation array 
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
            %    
            % Authors           
            % -------           
            %   Moses C. Nah    mosesnah@mit.edu
            % 
            %
            % Parameters
            % ----------
            %   (1) t_arr - The time array input
            %
            % Returns
            % -------          
            %   (1) act_sum - The array of sum of activation
            % ===========================================================================

            % Sum over the Basis Functions
            act_sum = sum( obj.calc_ith_arr( t_arr, 1:obj.N ), 1 );

        end

        function act_weighted = calc_whole_weighted_at_t( obj, t_arr, w_arr )
            % ===========================================================================
            % Descriptions
            % ------------
            %    fs = NonlinearForcingTerm.calc_whole_weighted_at_t( t_arr, w)arr )
            %         Calculate the whole activation of the N basis
            %         with the weighting array w_arr
            %    
            % Authors           
            % -------           
            %   Moses C. Nah    mosesnah@mit.edu
            % 
            % Parameters
            % ----------
            %   (1) t_arr - The time input, either scalar or array
            %
            %   (2) w_arr - weight array to be multiplied
            %               The array can be a multi-dimensional array
            %               n x N, where N is the number of basis functions.
            %                            n is the number of degrees of freedom
            %
            % Returns
            % -------
            %   (1) act_weighted - sum_{i=1}^{N}w_i phi_i( s( t ) )/sum_{i=1}^{N}phi_i( s( t ) )
            %       
            % ===========================================================================

            % t must be a row vector
            assert( isrow( t_arr ) )
            
            % The # of column of the weight array must be identical 
            % to the number of basis functions
            [ n, nc ] = size( w_arr );
            assert( nc == obj.N );
            
            % Initialize the nonlinear forcing term 
            act_weighted = zeros( n, length( t_arr ) );

            for i = 1 : n
                act_weighted( i, : ) = w_arr( i, : )' .* obj.calc_ith_arr( t_arr, 1:obj.N )/obj.calc_whole_at_t( t_arr );
            end
        end   

        function force_arr = calc_forcing_term( obj, t_arr, w_arr, t0i )
            % ===========================================================================
            % Descriptions
            % ------------
            %    fs = NonlinearForcingTerm.calc_forcing_term( t_arr, w_arr, t0i )
            %         Calculate the nonlinear forcing term
            %         with a time offset of t0i
            % ===========================================================================
    
            % First, we only need to calculate the forcing term between 
            % t0i and t0i + D, where D is the duration of the movement
            % D is simply the tau of canonical system. 
    
            % First, get the length of t_arr and w_arr, 
            % t must be row or column vector
            assert( isrow( t_arr ) )        
            Nt = length( t_arr );
    
            % Number of column of weight array must be equal to the number
            % of basis functions
            [ n, nc ] = size( w_arr );
            assert( nc == obj.N );
    
            % The time array's maximum value should be bigger than t0i
            assert( t0i <= max( t_arr ) );
    
            % Initialize the forcing term array
            force_arr = zeros( n, Nt );
    
            % Get the index for the time array
            % That is between t0i and t0i + D
            idx_arr = ( t_arr >= t0i && t_arr <= t0i + obj.cs.tau );
    
            % Calculate the forcing term 
            t_tmp = t_arr( idx_arr ) - t0i;
            force_arr( :, idx_arr ) = obj.cs.calc( t_tmp ) .* obj.calc_whole_weighted_at_t( t_tmp, w_arr );

        end

    end
end