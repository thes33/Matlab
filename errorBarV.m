function errorBarV(X, meanY, errorT, width, linewidth, color)
%function errorBarV(X, meanY, errorT, width, linewidth, color)
%
%  errorBarV:  Plot better vertical error bars, best used in cojunction with bar
%                 or plot functions.
%     INPUTS:
%        X = x coordinates for error bars (single column/row)
%        meanY - mean Y value around which bar is centered  (single column/row)
%        errorT - +/- this amount for error term, or [lower upper] boundary columns.
%        width (optional) - Desired width of the error bar in X units [Default 1].
%        linewidth (optional) - Line width parameter passed to line [Default 1].
%        color (optional) - Color parameter passed to line  [Default 'k'].
%
%  @Author: Eugene Brandewie 1/30/2013
%
%   See also:  PLOT, BAR, LINE
%


% Inputs
%---------------------------
if nargin < 6
    color = 'k';
end
if nargin < 5
    linewidth = 1;
end
if nargin < 4
    width = 1;
end
if nargin < 3
    error('Not enough input arguments: need mean and error term.');
end

if ((size(errorT,2)) > 2)
    error('Too many columns in error term: Lower and Upper bounds only.');
end

% Draw
%---------------------------
hw = width/2;

for xx = 1:length(X)
    
    if ((size(errorT,2)) == 2)
        lowerE = errorT(xx,1);
        upperE = errorT(xx,2);
    else
        lowerE = meanY(xx)-errorT(xx);
        upperE = meanY(xx)+errorT(xx);
    end
    
    line([X(xx) X(xx)], [lowerE upperE],'Linewidth',linewidth,'Color',color);
    line([X(xx)-hw X(xx)+hw], [lowerE lowerE],'Linewidth',linewidth,'Color',color);
    line([X(xx)-hw X(xx)+hw], [upperE upperE],'Linewidth',linewidth,'Color',color);
    
    
    
end
