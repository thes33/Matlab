function plotPolarInt(locations, values, varargin);
%function plotPolarInt(locations, values, varargin);
%
% Plots any range of values across azimuth locations in a colored
%   polar plot.
%
%  INPUTS:
%   locations - Array of evenly-spaced azimuth angles (from 0 to 360 degrees).
%   values(N,loc) - Variable values with rows representing condition types
%       and each column refers to a different location in the same order
%       as 'locations'.
%
% OPTIONS:
%   'hold' - Set to hold plots, useful for subplots.
%   'colormap', A - 64x3 array of color values to use for the color mapping.
%   'labels', C - Cell array of labels for each 'condition' or value type.
%   'showvalues' - Enables individual box labels to be printed.
%   'noboundaries' - Disables boundary lines between value type radii.
%
%
% @Author:  Eugene Brandewie (Mar 18th, 2019)





%% INPUT HANDLING
%======================================================
HOLD = 0;
COLOR_MAP = jet;
VALUE_LABELS = 0;
BOUNDARY_LINES = 1;

if nargin >= 2
    for ii = 1:nargin-2
        S = varargin{ii};
        if (~isnumeric(S))
            % PLOT
            if strcmpi(S,'hold')
                HOLD = 1;
            end
            % COLOR MAP
            if strcmpi(S,'colormap')
                COLOR_MAP = varargin{ii+1};
            end
            % LABELS
            if strcmpi(S,'labels')
                LABELS = varargin{ii+1};
                if (length(LABELS) ~= size(values,1)) 
                    error('Number of labels must match number of value types.');
                end
            end
            % VALUE LABELS
            if strcmpi(S,'showvalues')
                VALUE_LABELS = 1;
            end
            % BOUNDARY LINES
            if strcmpi(S,'noboundaries')
                BOUNDARY_LINES = 0;
            end
        end
    end
end



%% PREPARE INPUTS
%======================================================

degSep = 360/length(locations);
halfDegSep = degSep/2;

numCols = size(values,1);
colSep = 0.5;

% Default labels
if (~exist('LABELS','var'))
    LABELS = {};
    for (ll = 1:numCols)
        LABELS(ll) = {num2str(ll)};
    end
end


%% DRAW FIGURE
%======================================================

% Create Figure
if (~HOLD)
    figure();
end
hold on;
center = [0 0];
centerRad = (numCols/10);
radius = numCols + (centerRad) + colSep;

% Draw background
rectangle('Position',[-radius -radius (radius*2) (radius*2)],'Curvature',[1 1],'FaceColor',[1 1 1]);


% Plot outer circle
theta = linspace(5*pi/2, pi/2, 500)';
x = radius*cos(theta)+center(1);
y = radius*sin(theta)+center(2);
plot(x,y,'k','LineWidth',2);
axis equal; axis off;



%% PLOT DATA
%========================================================

% FOR each angle of arc
for (ll = 1:length(locations))
    loc = locations(ll);
    % FOR each value column
    for (vv = 1:numCols)
        % New Patch:
        X = []; Y = []; Z = [];
        [x,y] = pol2cart(deg2rad(loc + halfDegSep + 90), vv + centerRad - colSep);
        X = [X; x];
        Y = [Y; y];
        Z = [Z; values(vv,ll)];
        [x,y] = pol2cart(deg2rad(loc + halfDegSep + 90), vv + centerRad + colSep);
        X = [X; x];
        Y = [Y; y];
        Z = [Z; values(vv,ll)];
        
        [x,y] = pol2cart(deg2rad(loc - halfDegSep + 90), vv + centerRad + colSep);
        X = [X; x];
        Y = [Y; y];
        Z = [Z; values(vv,ll)];
        [x,y] = pol2cart(deg2rad(loc - halfDegSep + 90), vv + centerRad - colSep);
        X = [X; x];
        Y = [Y; y];
        Z = [Z; values(vv,ll)];
        patch(X,Y,Z,'EdgeColor','k');
        
        if (VALUE_LABELS)
            [X,Y] = pol2cart(deg2rad(loc+90),vv + centerRad);
            text(X, Y, num2str(values(vv,ll)), 'FontSize',10,'Color','k');
        end
    end
end



%% PLOT LABELS
%========================================================

% Plot Coordinate lines
line([0 0],[0 radius],'Color','k','LineStyle',':');
line([0 radius],[0 0],'Color','k','LineStyle',':');
line([0 0],[0 -radius],'Color','k','LineStyle',':');
line([0 -radius],[0 0],'Color','k','LineStyle',':');

[X,Y] = pol2cart(pi/4,radius);
line([0 X],[0 Y],'Color','k','LineStyle',':');
[X,Y] = pol2cart(-pi/4,radius);
line([0 X],[0 Y],'Color','k','LineStyle',':');
[X,Y] = pol2cart(3*pi/4,radius);
line([0 X],[0 Y],'Color','k','LineStyle',':');
[X,Y] = pol2cart(-3*pi/4,radius);
line([0 X],[0 Y],'Color','k','LineStyle',':');



% Plot inner 'head' circle
rectangle('Position',[-centerRad-colSep -centerRad-colSep (centerRad+colSep)*2 (centerRad+colSep)*2],'Curvature',[1 1],'FaceColor',[0.5 0.5 0.5]);



if (BOUNDARY_LINES)
    % Plot value type boundary circles
    %--------------------------------------
    % FOR each value column
    for (vv = 1:numCols)
        % Plot middle circles
        r = vv + centerRad + colSep;
        t = linspace(5*pi/2, pi/2, 500)';
        x = r*cos(t)+center(1);
        y = r*sin(t)+center(2);
        plot(x, y, 'k','LineWidth',2);
    end
end


% Labels
N = [10:10:180 -170:10:-10];
theta = [deg2rad(0) fliplr(deg2rad(N))];
N = [0 N];
r = 1.05*radius;
for (ii = 1:length(N))
    [X,Y] = pol2cart(theta(ii)+pi/2,r);
    str = strcat(num2str(N(ii),'%d'));
    text(X-(0.03*radius), Y, str, 'FontWeight','Bold');
end

% Value Type Labels
%------------------------------
% FOR each value column
for (vv = 1:numCols)
    [X,Y] = pol2cart(deg2rad(-88-halfDegSep),vv + centerRad + (colSep*.5));
    text(X, Y, LABELS(vv), 'FontSize',12,'FontWeight','Bold','Color','k');
end

% Side Color Bar
colorbar;
colormap(COLOR_MAP);








