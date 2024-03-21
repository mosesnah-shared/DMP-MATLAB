function activation = sigmoid_activation( t_arr, t_start, t_end, is_flip )
    % sigmoid_activation: A function which smoothly change from 0 to 1
    % 
    % Parameters:
    %   t_arr   - Time array
    %   t_start - The time at which the transition from 0 begins
    %   t_end   - The time at which the transition to 1 is completed
    %   is_flip - flip to change from 1 to 0
    %
    % Returns:
    %   activation - An array of the same size as t_arr with the smooth transition
    
    assert( isrow( t_arr) || iscolumn( t_arr) )

    % Ensure t_arr is a row vector for consistent operations
    if iscolumn( t_arr )
        t_arr = t_arr';
    end

    % Preallocate the return array with zeros
    activation = zeros( size( t_arr ) );
     
    % Indices where the transition occurs
    idx_arr = t_arr >= t_start & t_arr <= t_end;
    
    % Calculate the phase of the cosine function for the transition
    % The phase ranges from 0 to pi across the transition interval
    w = pi / (t_end - t_start);
    phi0 = t_start;

    % Apply the cosine function for the transition, mapping [1, -1] to [1, 0]
    activation( idx_arr ) = 0.5 * ( 1 + sin( w*( t_arr( idx_arr ) - phi0 ) - pi/2 ) );
    
    % Set the values to 0 after t_end
    activation( t_arr > t_end) = 1;

    if is_flip
        activation = 1-activation;
    end
end