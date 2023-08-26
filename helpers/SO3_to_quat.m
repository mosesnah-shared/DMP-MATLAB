function quat = SO3_to_quat( R )
% Method from the following reference:
% [REF] Allmendinger, Felix. Computational methods for the kinematic analysis of diarthrodial joints. 
%       Diss. Dissertation, Aachen, Techn. Hochsch., 2015, 2015.
%       Page 127, Algorithm 2

% A Quick check whether R is a rotation matrix 
assert( norm( R * R' - eye( 3 ), 'fro' ) <= 1e-9 );

% Calculate the Quaternion Matrix
R00 = trace( R );

% If trace is -1, then we need to do exception handling.

[~, k ] = max( abs( [ R00, R( 1, 1 ), R( 2, 2 ), R( 3, 3 ) ] ) );

if k == 1
    ek = 1/2 * sqrt( 1 + R00 );
else
    ek = 1/2 * sqrt( 1 + 2 * R( k-1, k-1 ) - R00 );
end


switch k-1
    case 0
        e0 = ek;
        e1 = 1/( 4 * ek ) * ( R( 3, 2 ) - R( 2, 3 ) );
        e2 = 1/( 4 * ek ) * ( R( 1, 3 ) - R( 3, 1 ) );
        e3 = 1/( 4 * ek ) * ( R( 2, 1 ) - R( 1, 2 ) );
    case 1
        e0 = 1/( 4 * ek ) * ( R( 3, 2 ) - R( 2, 3 ) );
        e1 = ek;
        e2 = 1/( 4 * ek ) * ( R( 2, 1 ) + R( 1, 2 ) );
        e3 = 1/( 4 * ek ) * ( R( 1, 3 ) + R( 3, 1 ) );        
    case 2
        e0 = 1/( 4 * ek ) * ( R( 1, 3 ) - R( 3, 1 ) );
        e1 = 1/( 4 * ek ) * ( R( 2, 1 ) + R( 1, 2 ) );
        e2 = ek;
        e3 = 1/( 4 * ek ) * ( R( 3, 2 ) + R( 2, 3 ) );        
    case 3        
        e0 = 1/( 4 * ek ) * ( R( 2, 1 ) - R( 1, 2 ) );
        e1 = 1/( 4 * ek ) * ( R( 1, 3 ) + R( 3, 1 ) );
        e2 = 1/( 4 * ek ) * ( R( 3, 2 ) + R( 2, 3 ) );               
        e3 = ek;
    otherwise
        % Will not happen just include
        assert( true )
end

   quat = [ e0, e1, e2, e3 ];
        
    if e0 < 0
       quat = -quat ;
    end

end

