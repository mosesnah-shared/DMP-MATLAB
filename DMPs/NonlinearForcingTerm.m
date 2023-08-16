classdef NonlinearForcingTerm < handle
    
    properties
        % Canonical system
        cs
        type
        
        % Number of Basis Function
        N

        % The center location and width of discrete and rhythmic movements 
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
            
            % Type of the nonlinear forcing term is naturally defined using the Canonical System
            obj.type = obj.cs.type;
            
            % Number of Basis Functions
            obj.N = N;

            % Setting the width and center locations
            % [2023.08.15] Requires improvement for the rhythmic movement case.
            % Case for Discrete movement followed the information from:
            % Saveriano, Matteo, et al. "Dynamic movement primitives in robotics: A tutorial survey." arXiv preprint arXiv:2102.03861 (2021).
            obj.c_arr = zeros( 1, N );
            obj.h_arr = zeros( 1, N );
            
            if obj.type == 0
                obj.c_arr = exp( -obj.cs.alpha_s/( obj.N - 1 ) * ( 0:(obj.N-1) ) );
                obj.h_arr( 1:end-1 ) = 1.0 ./ diff( obj.c_arr ).^2;
                obj.h_arr( end ) = obj.h_arr( end-1 );
            else
                obj.c_arr = 2*pi/obj.N * ( 0.5:1:obj.N );
                obj.h_arr = obj.N;
            end

        end

        
        function act = calc_ith( obj, t, i )
            
            % Number should be within the number of Basis Functions
            assert( i >= 1 && i <= obj.N )
            
            ci = obj.c_arr( i );
            hi = obj.h_arr( i );
            
            % The calculation of the canonical system
            s  = obj.cs.calc( t );
            
            % Calculate the activation of the i-th function
            % For discrete movement            
            if obj.type == 0
                act = exp( -hi * ( s - ci ).^2 );
                
            % For rhythmic movement                
            else
                act = exp( hi * ( cos( s - ci ) - 1 ) );
            end
        end
        
    end
end