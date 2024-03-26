function fig_save( hFig, figName )
% myColor Calling the list of colors for plot
%
% =============================================================== %
% [INPUT] 
%   (1) hFig 
%       -  The handle of the figure object that we are aimed to save.
%          The figure will be saved as a pdf file.
%
%   (2) figName:
%       -  The name of the figure pdf file.
%
% =============================================================== %
%
% =============================================================== %
% [OUTPUT] 
%   (-) ---
%       - FILL-IN
%
% =============================================================== %
%
% [REMARKS]  ADD DETAILS
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
% [   DATE   ]: 07-June-2020
% =============================================================== %

% =============================================================== %
% [DESCRIPTION]
%   Function for plotting a static 3D Plot 
%
% =============================================================== %

set( hFig, 'Units' , 'Inches' );
pos = get( hFig, 'Position' );
set( hFig, 'PaperPositionMode',   'Auto', ...
                  'PaperUnits', 'Inches', ...
                   'PaperSize', [ pos( 3 ), pos( 4 ) ] )
               
% saveas( hFig, figName + ".pdf" )
exportgraphics( hFig, [ figName, '.pdf' ],'ContentType','vector')


end