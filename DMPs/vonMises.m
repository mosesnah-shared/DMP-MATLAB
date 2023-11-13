function act = vonMises( s_arr, ci, hi )
% ===========================================================================
% Descriptions
% ------------
%    Von-Mises Function
%    f(s) = exp[ hi{ cos( s_arr - ci ) - 1 } ]
%
% Parameters
% ----------
%   (1) s_arr - A row array 
%   
%   (2) ci - Center of the von Mises Function
%
%   (3) hi - Width  of the von Mises Function
% 
% Returns
% -------
%   (1) act = exp[ hi{ cos( s_arr - ci ) - 1 } ]
%
% ===========================================================================

% t_arr must be a row vector
assert( isrow( s_arr ) )

% Center location and height of the von Mises Function
assert( isscalar( ci ) && isscalar( hi ) );

act = exp( hi * ( cos( s_arr - ci ) - 1 ) );

end

