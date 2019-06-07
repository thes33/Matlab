function reflectogram(theta, time, fs, varargin)
% function reflectogram(theta, time, fs, parameters...)
%
% Plots the reflections given by 'genBRIR' on a polar plot with timing
%   along the radius and azimuth along the circle.
%
%
% INPUTS:
%   theta = List of reflection azimuth angles (deg), outputted by 'genBRIR'.
%   time = List of reflection timings (samples), outputted by 'genBRIR'.
%       Assumed to be sorted by first to last in time.
%   fs = sampling frequency (Default: 48000)
%
% Optional parameters:
%      'hold' - Set to hold plots, useful for subplots.
%      'max',N - Maximum number of reflections to plot. [Default: 999]. -1 is all.
%      'single' - Mode to plot each individual reflection one at a time.
%           Use < and > to move between each reflection presentation.
%
%
%   REQUIRES SIGNAL PROCESSING TOOLBOX
%
% @Author: Eugene Brandewie - Apr 19, 2018


%% INPUT HANDLING
%======================================================
if (nargin < 2)
    error('Requires reflection azimuth and timing information.');
elseif (nargin < 3)
    fs = 48000;
end

% Parameter defaults:
HOLD = 0;
MAX = 1000;
SINGLE = 0;

% Optional settings:
if nargin >= 4
    for ii = 1:nargin-3
        S = varargin{ii};
        if (~isnumeric(S))
            % HOLD
            if strcmpi(S,'hold')
                HOLD = 1;
            end
            % MAX
            if strcmpi(S,'max')
                MAX = varargin{ii+1} + 1;
            end
            % SINGLE
            if strcmpi(S,'single')
                SINGLE = 1;
            end
        end
    end
end

% Show all reflections option
if (MAX == -1), MAX = length(theta); end
% Limit to total number inputed
if (MAX > length(theta)), MAX = length(theta); end

% Sort by timing
[timeMS,I] = sort(time);
thetaM = theta(I);
% Reduce set to MAX
timeMS = (timeMS(1:MAX)/fs).*1000;
thetaM = thetaM(1:MAX);




%% DRAW FULL FIGURE
%======================================================
if (~SINGLE)
    
    % Create Figure
    if (~HOLD)
        figure('Name','Reflectogram');
    end
    % Set figure parameters
    data = guidata(gcf);
    data.Reflectogram_currentReflection = length(thetaM);
    data.Reflectogram_maxReflections = length(thetaM);
    data.Reflectogram_SINGLE = false;
    data.Reflectogram_timeMS = timeMS;
    data.Reflectogram_thetaM = thetaM;
    guidata(gcf,data);
    
    % Render graph
    Reflectogram_RenderFrame(gcf);
    
else % ELSE SINGLE MODE
    %% SINGLE REFLECTION MODE
    %=====================================================
    
    % Create Figure
    FIG = figure('Name','Reflectogram', 'KeyPressFcn',@Reflectogram_KeyPressFcn);
    hold on;
    
    % Set figure parameters
    data = guidata(FIG);
    data.Reflectogram_currentReflection = 1;
    data.Reflectogram_maxReflections = length(thetaM);
    data.Reflectogram_SINGLE = true;
    data.Reflectogram_timeMS = timeMS;
    data.Reflectogram_thetaM = thetaM;
    guidata(FIG,data);
    
    % Render graph
    Reflectogram_RenderFrame(FIG);
    
end % END IF NOT SINGLE


end % END main function



%% REFLECTION KEYPRESS FUNCTION
%=============================================
function Reflectogram_KeyPressFcn(object, eventdata)
data = guidata(object);

% Previous: '<' or ',' Key
if (strcmp(eventdata.Key,'leftarrow') || strcmp(eventdata.Key,'comma'))
    if (data.Reflectogram_currentReflection == 1)
        data.Reflectogram_currentReflection = data.Reflectogram_maxReflections;
    else
        data.Reflectogram_currentReflection = max(1,data.Reflectogram_currentReflection-1);
    end
end

% Next: '>' or '.' Key
if (strcmp(eventdata.Key,'rightarrow') || strcmp(eventdata.Key,'period'))
    if (data.Reflectogram_currentReflection == data.Reflectogram_maxReflections)
        data.Reflectogram_currentReflection = 1;
    else
        data.Reflectogram_currentReflection = min(data.Reflectogram_maxReflections,data.Reflectogram_currentReflection+1);
    end
end

% Update data
guidata(object,data);

% Redraw frame
Reflectogram_RenderFrame(object);

end % END keypress function




%% DRAW REFLECTION FRAME
%=============================================
function Reflectogram_RenderFrame(object)
data = guidata(object);

timeMS = data.Reflectogram_timeMS;
thetaM = data.Reflectogram_thetaM;

% Clear figure if SINGLE
if (data.Reflectogram_SINGLE)
    clf;
end
hold on;

% Set both circle radius'
center = [0 0];
radius = max(timeMS);
centerRad = min(timeMS)-1;
%centerRad = radius * 0.1;

% Draw background
rectangle('Position',[-radius -radius (radius*2) (radius*2)],'Curvature',[1 1],'FaceColor',[1 1 1]);


% Plot outer circle
T = linspace(5*pi/2, pi/2, 500)';
x = radius*cos(T)+center(1);
y = radius*sin(T)+center(2);
plot(x,y,'k','LineWidth',2);
axis equal; axis off;



%% PLOT DATA
%========================================================

% Plot direct waveform
[x,y] = pol2cart(deg2rad(-thetaM(1) + 90), timeMS(1));
plot(x,y,'ok','MarkerSize',10,'MarkerFaceColor','r');

% FOR each reflection
if (data.Reflectogram_currentReflection > 1)
    for (ll = 2:data.Reflectogram_currentReflection)
        loc = thetaM(ll);
        tim = timeMS(ll);
        
        [x,y] = pol2cart(deg2rad(-loc + 90), tim);
        % IF SINGLE, highlight current reflection
        if (data.Reflectogram_SINGLE && ...
                ll == data.Reflectogram_currentReflection)
            plot(x,y,'ok','MarkerSize',10,'MarkerFaceColor','k');
        else
            plot(x,y,'ok','MarkerSize',4,'MarkerFaceColor','k');
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


% Plot middle circles
rr = [0.2*radius 0.4*radius 0.6*radius 0.8*radius];
for ii = 1:length(rr)
    r = rr(ii);
    if (rr(ii) > centerRad)
        x = r*cos(T)+center(1);
        y = r*sin(T)+center(2);
        plot(x, y, 'r:');
    end
end


% Plot inner 'head' circle
rectangle('Position',[-centerRad-1 -centerRad-1 (centerRad+1)*2 (centerRad+1)*2],'Curvature',[1 1],'FaceColor',[0.5 0.5 0.5]);


% Outer Azimuth Labels
N = [10:10:180 -170:10:-10];
T = [deg2rad(0) fliplr(deg2rad(N))];
N = [0 N];
r = 1.05*radius;
for (ii = 1:length(N))
    [X,Y] = pol2cart(T(ii)+pi/2,r);
    str = strcat(num2str(N(ii),'%d'));
    text(X-(0.03*radius), Y, str, 'FontWeight','Bold');
end

% Radius Timing Labels
rr = [0.2 0.4 0.6 0.8 1.0];
for ii = 1:length(rr)
    if (rr(ii)*radius > centerRad)
        [X,Y] = pol2cart(deg2rad(45),rr(ii)*radius);
        str = strcat(num2str(round(rr(ii)*max(timeMS)),'%d'));
        text(X, Y, str, 'FontWeight','Bold','Color',[1 0.3 0.3]);
    end
end

% IF SINGLE, show current reflection number  and time
if (data.Reflectogram_SINGLE)
    % Reflection Number
    [X,Y] = pol2cart(deg2rad(-40)+(pi/2),radius*1.3);
    text(X,Y,num2str(data.Reflectogram_currentReflection-1));
    % Timing
    [X,Y] = pol2cart(deg2rad(-42)+(pi/2),radius*1.2);
    text(X,Y,sprintf('%4.2f',data.Reflectogram_timeMS(data.Reflectogram_currentReflection)));
end



guidata(object,data);
end % END render frame function




