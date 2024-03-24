function dq_traj = data_diff( q_traj, t_arr )

% Data must be 2D 
assert( ismatrix( q_traj ) );

% Time array 
assert( iscolumn( t_arr) || isrow( t_arr ) );

% Assert that the data is appended along the column
[ nr, nc ] = size( q_traj );
assert( nc > nr );
assert( nc == length( t_arr ) );

% Differentiate and append the final point
dq_traj = diff( q_traj, 1, 2 )./ diff( t_arr );
dq_traj( :, end + 1 ) = dq_traj( :, end );

end