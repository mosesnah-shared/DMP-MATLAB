classdef CanonicalSystem < handle
    
    properties       
        
        % Type of the canonical system, either discrete (0), rhythmic (1)
        type 
        
        % The time constant of the canonical system
        % For Discrete movement: Duration of the movement
        % For Rhythmic movement: Period of the movement divided by 2pi
        tau 
        
        % Must be Positive value
        alpha_s
    end

    methods
        function obj = CanonicalSystem( type, tau, alpha_s )
            % ===========================================================================
            % Descriptions
            % ------------
            %    Default Constructor of the Canonical System
            %    cs = CanonicalSystem( type, tau, alpha_s )
            %
            % Author           
            % -------           
            %   Moses C. Nah    mosesnah@mit.edu
            % 
            %
            % Parameters
            % ----------
            %   (1) type - 'discrete' or 'rhythmic'
            % 
            %   (2) tau - The time constant of the canonical system
            %             For 'discrete' movement: Duration of the movement
            %             For 'rhythmic' movement: Period of the movement divide by 2pi
            %
            %   (3) alpha_s - positive constant, if 'rhythmic', then 
            %                 value is ignored
            %
            % ===========================================================================
            
            % Type input should be either 'discrete' or 'rhythmic'
            type = lower( type );
            assert( strcmp( type, "discrete" ) || strcmp( type, "rhythmic" ) )
            
            % Discrete (0) or Rhythmic (1)
            if strcmp( type, "discrete" )
                obj.type = 0;
            else
                obj.type = 1;
            end
            
            % Tau and alpha_s should be positive values
            assert( tau > 0 && alpha_s > 0 );
            obj.tau     = tau; 
            obj.alpha_s = alpha_s;

        end

        function s = calc( obj, t )
            % ===========================================================================
            % Descriptions
            % ------------
            %    Calculating the Canonical System
            %
            % Parameters
            % ----------
            %   (1) t - time (sec)
            %           accepts array input
            % 
            % Returns
            % -------
            %   (1) s - the calculataion of s(t)
            %           If discrete (0): s(t) = exp( -alpha_s/tau t )
            %           If rhythmic (1): s(t) = mod( t/tau, 2pi )
            %
            % ===========================================================================
            
            % If Discrete
            if obj.type == 0 
                s = exp( -obj.alpha_s/obj.tau * t );

            % If Rhythmic
            elseif obj.type == 1
                s = mod( t/obj.tau, 2*pi );

            % IF not, then should be halted
            else
                error( '[Wrong input] type should be either 0 or 1 but %d is defined as type', obj.type )
            end
            
        end
    end
end