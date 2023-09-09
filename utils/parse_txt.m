function output = parse_txt( dir, n_skip )

% Open the file
f = fopen( dir );

% Save the data as an output matrix output 
tline = fgetl( f );

% Number of lines to skip, should be integer
assert( n_skip >= 0 && rem( n_skip, 1 ) == 0 );

i = 1;
while ischar( tline )

    % Skip the lines
    if i < ( n_skip + 1 )
        i = i + 1;
        tline = fgetl( f );    
        continue
    end

    tmp = str2double( regexp( tline ,'[-+]?\d+\.?\d*','match') );
    output( i, : ) = tmp;
    i = i + 1;
    tline = fgetl( f );
end

end

