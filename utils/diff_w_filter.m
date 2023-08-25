function dq_traj = diff_w_filter( q_traj, filter_name, N )
% Differential q_traj with the filter, 

assert( filter_name == "gaussian" || filter_name == "none" );
assert( N>= 1 )
if filter_name == "none":
    is_filter = false;
else
    is_filter = true;   
end

% Get the size of trajectory
% The trajectory should be long along the columns, i.e., nq x N array
% The column should be longer than the row 
[ nr, nc ] = size( q_traj );

assert( nc > nr );

dq_traj = diff( q_traj, 1, 2 );

% Append the final value
dq_traj( :, end + 1 ) = dq_traj( :, end );

for i = 1 : nr
    dq_traj( i, : ) = smoothdata( dq_traj( i, : ), filter_name, N );
end

end

end