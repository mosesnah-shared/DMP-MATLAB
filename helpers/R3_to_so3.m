function so3 = R3_to_so3( vec )

assert( iscolumn( vec ) || isrow( vec ) );
assert( length( vec ) == 3 );

so3 = [       0, -vec( 3 ),  vec( 2 ), ...
         vec( 3 ),        0, -vec( 1 ), ... 
        -vec( 2 ), vec( 1 ),         0 ];
         
end