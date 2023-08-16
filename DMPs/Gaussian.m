function act = Gaussian( s, ci, hi )
% ===========================================================================
% Descriptions
% ------------
%    Non-normalized Gaussian Function
%    f(s) = exp( -hi(s-ci)^2 )
%
% Parameters
% ----------
%   (1) s  - value of x, either could be an array or scalar
%   
%   (2) ci - Center of the Gaussian
%
%   (3) hi - Width of the Gaussian
% 
% Returns
% -------
%   (1) act = exp( -hi(s-ci)^2 )
%
% ===========================================================================

assert( isscalar( ci ) && isscalar( hi ) );

act = exp( -hi * ( s - ci ).^2 );

end

