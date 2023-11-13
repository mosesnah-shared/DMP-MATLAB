function act = Gaussian( s_arr, ci, hi )
% ===========================================================================
% Descriptions
% ------------
%    Non-normalized Gaussian Function
%    f(s) = exp( -hi(s-ci)^2 )
%
% Parameters
% ----------
%   (1) s_arr - A row array 
%   
%   (2) ci    - Center of the Gaussian
%
%   (3) hi    -  Width of the Gaussian
% 
% Returns
% -------
%   (1) act = exp{ -hi ( s - ci )^2 }
%
% ===========================================================================

% t_arr must be a row vector
assert( isrow( s_arr ) )

% Center location and height of the Gaussian Function
assert( isscalar( ci ) && isscalar( hi ) );

act = exp( -hi * ( s_arr - ci ).^2 );

end

