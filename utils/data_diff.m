function dq_traj = data_diff( q_traj )

% Data must be 2D 
assert( length( size( q_traj ) ) == 2 );

% Assert that the data is appended along the column
[ nr, nc ] = size( q_traj );
assert( nc > nr );

% Differentiate and append the final point
dq_traj = diff( q_traj, 1, 2 );
dq_traj( :, end + 1 ) = dq_traj( :, end );

end