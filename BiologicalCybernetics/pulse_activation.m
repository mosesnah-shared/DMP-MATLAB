function activation = pulse_activation( t_arr, t_start, t_end )
% sigmoid_activation: A function which smoothly change from 0 to 1 and 1 to 0 
% 
% Parameters:
%   t_arr   - Time array
%   t_start - The time at which the transition from 0 begins
%   t_end   - The time at which the transition to 0 is completed
%   is_flip - flip to change from 1 to 0
%
% Returns:
%   activation - An array of the same size as t_arr with the smooth transition

% Initialize the transition array with zeros
assert( isrow( t_arr ) || iscolumn( t_arr ) );

if iscolumn( t_arr )
    t_arr = t_arr';
end

activation = zeros( size( t_arr ) );

% Determine the indices where the transition occurs
idx_arr = t_arr >= t_start & t_arr <= t_end;

% Calculate the phase of the sine wave for the transition
% The phase ranges from 0 to pi across the transition interval
phase = 2*pi * ( t_arr( idx_arr ) - t_start ) / ( t_end - t_start );

% Apply the sine function for the transition, scaling its output from 0 to 1
activation( idx_arr ) = 1-( 1 + cos( phase ) ) / 2;

end

