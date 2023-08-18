function output = parse_txt( dir )

% Open the file
f = fopen( dir );

% Save the data as an output matrix output 
tline = fgetl( f );

i = 1;
while ischar( tline )
    tmp = str2double( regexp( tline ,'[-+]?\d+\.?\d*','match') );
    output( i, : ) = tmp;
    i = i + 1;
    tline = fgetl( f );
end

end

