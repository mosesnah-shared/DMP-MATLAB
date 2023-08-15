classdef CanonicalSystem < handle
    
    properties
        type 
        tau 
        alpha_s
    end

    methods
        function obj = CanonicalSystem( type, tau, alpha_s )

            type = lower( type );
            assert( strcmp( type, "discrete" ) || strcmp( type, "rhythmic" ) )
            

            if strcmp( type, "discrete" )
                obj.type = 0;
            else
                obj.type = 1;
            end

            assert( tau > 0 && alpha_s > 0 );
            obj.tau = tau; 

            % alpha_s is only used for discrete movement
            % and will be ignored  for rhythmic movement.
            obj.alpha_s = alpha_s;

        end

        function s = calc( obj, t )
            if obj.type == 0
                s = exp( -obj.alpha_s/obj.tau * t );
            else
                s = mod( 1.0/obj.tau * t, 2*pi );
            end
        end
    end
end