function plotPolarSpectElev(locations, spect, freqs, varargin);
%function plotPolarSpectElev(locations, spect, freqs, varargin);
%
% Plots the frequency spectrum across elevation locations in a colored
%   polar plot.
%
%  INPUTS:
%   locations - Array of elevation angles (from 0 to 360).
%           0 is directly ahead, negative angles below, positive above
%   spect(f,l) - Spectrum energy with rows represented by 'freqs' bands
%       and each column refers to a different location in the same order 
%       as 'locations'.
%   freqs - Array of center frequencies for each band in spect.
%
%   'hold' - Set to hold plots, useful for subplots.
%
%  REQUIRES: fourierBands, SIGNAL PROCESSING TOOLBOX
%
% @Author:  Eugene Brandewie (Feb 28th, 2019)





%% INPUT HANDLING
%======================================================
HOLD = 0;

if nargin >= 4
    for ii = 1:nargin-3
        S = varargin{ii};
        if (~isnumeric(S))
            % PLOT
            if strcmpi(S,'hold')
                HOLD = 1;
            end
        end
    end
end



%% PREPARE INPUTS
%======================================================

% FOR each location angle
for ll = 1:length(locations)
    if (ll-1 < 1),  LB = 0;
    else LB = locations(ll) - locations(ll-1); end
    if (ll+1 > length(locations)), UB = 0; 
    else UB = locations(ll+1) - locations(ll); end    
    halfDegSep(ll) = max(abs(UB),abs(LB))/2;
    LowerBound(ll) = locations(ll) - halfDegSep(ll);
    UpperBound(ll) = locations(ll) + halfDegSep(ll);
end


%% DRAW FIGURE
%======================================================

% Create Figure
if (~HOLD)
    figure();
end
hold on;
center = [0 0];
centerRad = (length(freqs)/10);
radius = length(freqs)+(centerRad);

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
    
    
    % New Patch:
    X = []; Y = []; Z = [];
    % FOR each frequency band
    for (ff = 1:length(freqs))
        [x,y] = pol2cart(deg2rad(UpperBound(ll)), ff+centerRad);
        X = [X; x];
        Y = [Y; y];
        Z = [Z; spect(ff,ll)];
    end
    % FOR each frequency band in reverse
    for (ff = length(freqs):-1:1)
        [x,y] = pol2cart(deg2rad(LowerBound(ll)), ff+centerRad);
        X = [X; x];
        Y = [Y; y];
        Z = [Z; spect(ff,ll)];
    end
    patch(X,Y,Z,'EdgeColor','k');    
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


% Plot middle circles
rr = [0.2*radius 0.4*radius 0.6*radius 0.8*radius];
for ii = 1:length(rr)
    r = rr(ii);
    x = r*cos(theta)+center(1);
    y = r*sin(theta)+center(2);
    plot(x, y, 'r:');
end


% Plot inner 'head' circle
rectangle('Position',[-centerRad-1 -centerRad-1 (centerRad+1)*2 (centerRad+1)*2],'Curvature',[1 1],'FaceColor',[0.5 0.5 0.5]);


% Labels
N = [10:10:180 -170:10:-10];
theta = [deg2rad(0) (deg2rad(N))];
N = [0 N];
r = 1.05*radius;
for (ii = 1:length(N))
    [X,Y] = pol2cart(theta(ii),r);
    str = strcat(num2str(N(ii),'%d'));
    text(X-(0.03*radius), Y, str, 'FontWeight','Bold');
end

rr = [0.2 0.4 0.6 0.8 1.0];
for ii = 1:length(rr)
    [X,Y] = pol2cart(45,rr(ii)*radius);
    str = strcat(num2str(round(rr(ii)*max(freqs),-2),'%d'));
    text(X, Y, str, 'FontWeight','Bold','Color',[1 0.3 0.3]);
end
colorbar;









