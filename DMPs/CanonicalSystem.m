 classdef CanonicalSystem < handle
    
    properties       
        
        % Type of the canonical system, either discrete (0), rhythmic (1)
        type 
        
        % The time constant of the canonical system
        % For Discrete movement: Duration of the movement
        % For Rhythmic movement: Period of the movement divided by 2pi
        tau 
        
        % For Discrete movement: must be positive value
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
            % Parameters
            % ----------
            %   (1) type - 'discrete' or 'rhythmic'
            % 
            %   (2) tau - The time constant of the canonical system
            %             For 'discrete' movement: Duration of the movement
            %             For 'rhythmic' movement: Period of the movement divide by 2pi
            %
            %   (3) alpha_s - positive constant. 
            %                 If 'rhythmic', then value is ignored
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

        function s_arr = calc( obj, t_arr )
            % ===========================================================================
            % Descriptions
            % ------------
            %    Calculating the Canonical System value at t
            %
            % Parameters
            % ----------
            %   (1) t_arr - Time array (sec), as a row vector.
            %           
            % 
            % Returns
            % -------
            %   (1) s_arr - Calculation of s(t)
            %               If discrete (0): s(t) = exp( -alpha_s/tau t )
            %               If rhythmic (1): s(t) = mod( t/tau, 2pi )
            %
            % ===========================================================================
            
            % Time array t_arr should be a row vector
            assert( isrow( t_arr ) )

            % If Discrete
            if obj.type == 0 
                s_arr = exp( -obj.alpha_s/obj.tau * t_arr );

            % If Rhythmic
            elseif obj.type == 1
                s_arr = mod( t_arr/obj.tau, 2*pi );

            % IF not, then should be halted
            else
                error( '[Wrong input] type should be either 0 or 1 but %d is defined as type', obj.type )
            end
            
        end
    end
end