function func_addSubfolders( varargin )
% ===============================================================================
% func_addSubfolders - Adding the list provided by varargin to the path
%
% Authors                         Email
%   Johannes Lachner              jlachner@mit.edu
%   Moses C. Nah                  mosesnah@mit.edu
% 
% Input 
%   [1] varargin - the cell list of strings
%
% ===============================================================================

% Assert that the varargin should all be a cell string
assert( iscellstr( varargin ) )

N = length( varargin );

for i = 1: N
    
    cur_path = [ pwd, '/' , varargin{ i } ];
    addpath( cur_path );
    addpath( genpath( cur_path ) );

end
