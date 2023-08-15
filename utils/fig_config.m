function fig_config( varargin )
% myFigureConfig  Setting the default figure/axis configuration 
%
% =============================================================== %
% [INPUT] varagin: getting the whole list of inputs given
%        
% [PROPERTIES]                                                  [DEFAULT]
%    fontSize: the default size of the figure font.                10
%   lineWidth: the default width of the line.                       5
%  markerSize: the default size of the marker.                     20
% =============================================================== %
%
% =============================================================== %
% [OUTPUT]    None
%
%
% =============================================================== %
%
% [EXAMPLES] (1) mySetFigureConfig( 'fontSize', 10 )
%            (2) mySetFigureConfig( 'fontSize', 10, 'markerSize', 10 )
%
%
% =============================================================== %
%
% =============================================================== %
% SEE ALSO testHelpFunction 
%
% =============================================================== %
% 
% [CREATED BY]: Moses C. Nah
% [DATE]      : 07-June-2020
%
% =============================================================== %

% =============================================================== %
% [REFERENCES]
%
% [Primitive Objects Properties]: https://www.mathworks.com/help/matlab/graphics-object-properties.html
% [Figure Properties]           : https://www.mathworks.com/help/matlab/ref/matlab.ui.figure-properties.html
%
%
% =============================================================== %


p = inputParser( );                                                        % Case sensitive is OFF By default
                                                                           % [REF] https://de.mathworks.com/help/matlab/ref/inputparser.html

% [TIP] We can add a function handle for a simple input-type check
ckc = @(x) ( isnumeric(x) && ( x > 0 ) );

addParameter( p,      'fontSize',  10, ckc );
addParameter( p,     'lineWidth',  5 , ckc );
addParameter( p,    'markerSize',  20, ckc );   
addParameter( p, 'axesLineWidth', 1.8, ckc );

parse( p, varargin{ : } )

r = p.Results;

fs  = r.fontSize;   
lw  = r.lineWidth;
ms  = r.markerSize;
alw = r.axesLineWidth;

set( 0, 'defaultTextfontSize'               ,     1.6*fs        );
set( 0, 'defaultTextInterpreter'            ,     'latex'       );
set( 0, 'defaultLegendInterpreter'          ,     'latex'       );
set( 0, 'defaultLineLinewidth'              ,       lw          );
set( 0, 'defaultAxesLinewidth'              ,        alw        );
set( 0, 'defaultLineMarkerSize'             ,        ms         );
set( 0, 'defaultAxesTickLabelInterpreter'   ,     'latex'       );  
set( 0, 'defaultAxesfontSize'               ,        fs         );
set( 0, 'defaultAxesXGrid'                  ,       'on'        );
set( 0, 'defaultAxesYGrid'                  ,       'on'        );
set( 0, 'defaultAxesZGrid'                  ,       'on'        );
set( 0, 'defaultAxesBox'                    ,       'on'        );
set( 0, 'defaultFigureWindowStyle'          ,     'normal'      );
set( 0, 'defaultFigureUnits'                ,     'normalized'  );
set( 0, 'defaultFigurePosition'             ,     [0 0 1 1]     );
set( 0, 'defaultFigureColor'                ,     [1 1 1  ]     );

set( 0, 'defaultFigureCreateFcn'            , @( fig, ~ )addToolbarExplorationButtons( fig ) )
set( 0, 'defaultAxesCreateFcn'              , @(  ax, ~ )set( ax.Toolbar, 'Visible', 'off' ) )



end