function act = vonMises( s, ci, hi )
% ===========================================================================
% Descriptions
% ------------
%    Von-Mises Function
%    f(s) = exp[ hi{ cos( s - ci ) - 1 } ]
%
% Parameters
% ----------
%   (1) s  - value of x, either could be an array or scalar
%   
%   (2) ci - Center of the von Mises Function
%
%   (3) hi - Width  of the von Mises Function
% 
% Returns
% -------
%   (1) act = exp[ hi{ cos( s - ci ) - 1 } ]
%
% ===========================================================================

assert( isscalar( ci ) && isscalar( hi ) );

act = exp( hi * ( cos( s - ci ) - 1 ) );

end

