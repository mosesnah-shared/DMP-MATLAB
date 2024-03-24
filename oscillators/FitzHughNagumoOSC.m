classdef FitzHughNagumoOSC
    % A FitzHugh-Nagumo Oscillator
    % 2 State Variables
    % [REF] Nagumo, Jinichi, Suguru Arimoto, and Shuji Yoshizawa. 
    % "An active pulse transmission line simulating nerve axon." 
    % Proceedings of the IRE 50.10 (1962): 2061-2070.

    properties
        a      % System Property #1, Refer below
        b      % System Property #2, Refer below
        c      % System Property #3, Refer below
    end

    methods
        function obj = FitzHughNagumoOSC( a, b , c)
            
            assert( b > 0 && b < 1 );
            assert( c^2 > b  );
            assert( ( 1 > a ) && ( a > 1-2/3*b ) );

            obj.a = a;
            obj.b = b;
            obj.c = c;
         
        end

        function x_new = step( obj, dt, x_old, stim, input )
            assert( isrow( x_old ) || iscolumn( x_old ) );
            assert( isrow( input ) || iscolumn( input ) );
            assert( isscalar( stim )  );

            assert( length( x_old ) == 2 );
            assert( length( input ) == 2 );
            
            if isrow( x_old )
                x_old = x_old';
            end

            v = x_old( 1 ); w = x_old( 2 );

            dv =    obj.c * ( v + w - 1/3*v^3 + stim );
            dw = -1/obj.c * ( v - obj.a + obj.b*w );

            x_new = x_old + [ dv; dw]*dt + input;
        end
    end
end