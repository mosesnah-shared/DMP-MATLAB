classdef CanonicalSystem < handle
    
    properties       
        
        % Type of the canonical system, either discrete (0), rhythmic (1)
        type 
        
        % The time constant of the canonical system
        % For Discrete movement: Duration of the movement
        % For Rhythmic movement: Period of the movement divide by 2pi
        tau 
        
        % Positive value
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
            % Authors           
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
            %   (3) alpha_s - positive constant
            %
            % ===========================================================================
            
            % Type input should be either discrete or rhythmic
            type = lower( type );
            assert( strcmp( type, "discrete" ) || strcmp( type, "rhythmic" ) )
            
            % Discrete (0) or Rhythmic (1)
            if strcmp( type, "discrete" )
                obj.type = 0;
            else
                obj.type = 1;
            end
            
            % Tau and alpha_s 
            assert( tau > 0 && alpha_s > 0 );
            obj.tau = tau; 
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
            % 
            % Returns
            % -------
            %   (1) s - the calculataion of s(t)
            %           If discrete: s(t) = exp( -alpha_s/tau t )
            %           If rhythmic: s(t) = mod( t/tau, 2pi )
            %
            % ===========================================================================
            
            % Discrete
            if obj.type == 0 
                s = exp( -obj.alpha_s/obj.tau * t );
                
            % Rhythmic    
            else             
                s = mod( t/obj.tau, 2*pi );
            end
            
        end
    end
end