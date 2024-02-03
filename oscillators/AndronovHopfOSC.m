classdef AndronovHopfOSC < handle
    % An Andronov-Hopf Oscillator
    % 2 State Variables

    properties
        r      % Radius of the Limit Cycle
        w      % Angular Velocity          
    end

    methods
        function obj = AndronovHopfOSC( w, r )
            assert( isscalar( w ) && w > 0 );
            assert( isscalar( r ) );
            obj.w = w;
            obj.r = r;

        end

        function x_new = step( obj, dt, x_old, input )
            assert( isrow( x_old ) || iscolumn( x_old ) );
            assert( isrow( input ) || iscolumn( input ) );

            assert( length( x_old ) == 2 );
            assert( length( input ) == 2 );
            
            if isrow( x_old )
                x_old = x_old';
            end

            x1 = x_old( 1 ); x2 = x_old( 2 );
            A = zeros( 2, 2 );
            A( 1, 1 ) =  obj.r^2-x1^2-x2^2;
            A( 1, 2 ) = -obj.w;
            A( 2, 1 ) =  obj.w;
            A( 2, 2 ) =  obj.r^2-x1^2-x2^2;
          
            x_new = x_old + dt * A * x_old + input;
        end
    end
end