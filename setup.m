% [Project] Exp[licit] 
% - The robotic simulator based on differential geometry
% You must run this setup.m file to run the .m scripts.
%
% Authors                       Email
%   [1] Johannes Lachner        jlachner@mit.edu
%   [2] Moses C. Nah            mosesnah@mit.edu
%
%
% The code is heavily commented. A famous quote says:
% "Code is read more often than it is written"
%           - Guido Van Rossum, the Creator of Python

%% Cleaning up + Environment Setup
clear; close all; clc;

% Include all the subdirectories
add_folders( 'DMPs', 'examples', 'trajectory', 'utils', 'export_fig', ... 
             'GeometryLibrary/MATLAB', 'example-contraction', 'ThesisImages', ...
             'oscillators', 'synchronization' );

% Also run the setup script under Explicit-MATLAB
cd( 'Explicit-MATLAB' );
setup
cd( '..' );