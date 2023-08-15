classdef NonlinearForcingTerm

    properties
        % Canonical system
        cs
        
        % Number of Discrete Movements
        N

        % [Discrete] Goal and Initial Location
        g
        y0

        % [Rhythmic] Scaling Coefficient
        r

        % The center location and width of discrete and rhythmic movements 
        c_arr = zeros( 1, N );
        h_arr = zeros( 1, N );

    end

    methods
        function obj = NonlinearForcingTerm( cs, N, g, y0, r )
            % The canonical system of the Nonlinear Forcing term 
            obj.cs = cs; 
            
            % Number of Basis Functions
            obj.N = N;

            % [Discrete Movement] Initial and Goal Location 
            obj.g  = g;
            obj.y0 = y0;
            
            % [Rhythmic Movement] Scaling Factor
            obj.r  = r;

            % Setting the width and center locations
            if obj.cs.type == 0
                obj.c_arr = exp( -obj.cs.alpha_s/( obj.N - 1 ) * [ 0:obj.N-1 ] ) );
                obj.h_arr( 1:end-1 ) = 1 / diff( obj.c_arr ).^2;
                obj.h_arr( end ) = obj.h_arr( end-1 );
            else
                obj.c_arr = 2*pi/obj.N * [ 0.5:1:obj.N ];
                obj.h_arr = obj.N;
            end

        end


        function outputArg = calc_ith( obj, t )
            
        end
    end
end