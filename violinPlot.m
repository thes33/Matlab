function [STATS] = violinPlot(X, Y, varargin);
%function violinPlot(X, Y, parameters...);
%
%  violinPlot:  Plot violin plots with error bars for each column of data.
%
%     INPUTS:
%        X = x coordinates for each column (center of violins).
%        Y = Y data value columns to determine means and errors.
%     
%     OPTIONAL INPUTS:
%        'bins', N = Number of bins for frequency analysis [Default: 5].
%        'colors', C = Color of the underlying violin patch.  [Default 'k'].
%            Can have a value for each data column as cell array, Ex: {'k','b','r'}.
%        'width', N = Desired maximum width of the violin plot densities. [Default 0.8].
%        'points' = Plots individual data points as X's. Outliers in red.
%        'box' = Plots the interquartile range box (1.5 x IQR) for outliers.
%        'outliers' = Ignores outliers (includes them in statistics).
%            Outliers defined as outside 1.5 x interquartile range.
%
%     OPTIONAL INPUTS:
%
%        STATS = Structure of statistical properties for each data column.
%             .MEAN = Arithmetic mean of each data column.
%             .MED = Median of each data column.
%             .IQR = Interquartile range of each data column.
%             .Q1 = First quartile median of each data column.
%             .Q3 = Third quartile median of each data column.
%             .LQ = Q1 - 1.5 * IQR, lower bound for outliers.
%             .UQ = Q3 + 1.5 * IQR, upper bound for outliers.
%             .LB = 95% confidence interval lowerbounds.
%             .UB = 95% confidence interval upperbounds.
%
%  REQUIRES: confidenceInterval
%
%  @Author: Eugene Brandewie 12/07/2017
%
%  VERSION: 
%   1.0 - Initial Release
%   1.1 - Added 'box' option and statistical output.
%   1.2 - Removed reliance on nan suite. Replaced with 'omitnan' tags.
%


%% PROCESS INPUTS
%===================================
% Set defaults
color = 'k';
width = 0.8;
bins = 5;
IND_POINTS = false;
IGNORE_OUTLIERS = false;
PLOT_BOX = false;

% Error control
if nargin < 2
    error('Not enough input arguments: need data columns.');
end
if (length(X) ~= size(Y,2))
    error('Size of X must equal number of data columns (Y).');
end

% FOR each parameter
if nargin > 2
    for ii = 1:(nargin-2)
        Q = varargin{ii};
        if (~isnumeric(Q))
            % Color
            if (strcmp(Q,'colors'))
                color = varargin{ii+1};                
            end
            % Width
            if (strcmp(Q,'width'))
                width = varargin{ii+1};                
            end
            % Bins
            if (strcmp(Q,'bins'))
                bins = varargin{ii+1};                
            end
            % Individual Data Points
            if (strcmp(Q,'points'))
                IND_POINTS = true;        
            end
            % Ignore Outliers
            if (strcmp(Q,'outliers'))
                IGNORE_OUTLIERS = true;        
            end
            % Box
            if (strcmp(Q,'box'))
                PLOT_BOX = true;        
            end
        end
    end
end

if (bins < 1)
    error('Number of bins must be at least 1.');
end
% Limit bin number to sample size
if (bins > size(Y,1))
    warning(['Bins reduced to sample size (' num2str(size(Y,1)) ')']);
end
% Ensure color values for each data column
if (length(color) < size(Y,2) || length(color) < 2)
    C = color;
    color = {};
    for (ii=1:size(Y,2))
        color(ii) = { C };
    end
end


%% CALCULATE STATISTICS
%=====================================

% FOR each column of data
for (dd = 1:size(Y,2))
    
    % Process Outliers
    if (~IGNORE_OUTLIERS)
        % Determine Outliers
        Med = median(Y(:,dd),'omitnan');
        Q1(dd) = median(Y(find(Y(:,dd)<Med),dd),'omitnan');
        Q3(dd) = median(Y(find(Y(:,dd)>Med),dd),'omitnan');
        
        IQR(dd) = 1.5 * (Q3(dd) - Q1(dd));
        L(dd) = Q1(dd) - IQR(dd);
        U(dd) = Q3(dd) + IQR(dd);
        
        outliers(dd) = {find(Y(:,dd) > U(dd) | Y(:,dd) < L(dd))};
        normdata(dd) = {find(Y(:,dd) < U(dd) & Y(:,dd) > L(dd))};
        if (~isempty(outliers{dd}))
            warning([num2str(length(outliers{dd})) ' outliers found in column ' num2str(dd) '.']);
        end
    else % Else ignore outliers
        Med = median(Y(:,dd),'omitnan');
        normdata(dd) = {1:length(Y(:,dd))};
        outliers(dd) = {[]};
        
        Q1(dd) = median(Y(find(Y(:,dd)<Med),dd),'omitnan');
        Q3(dd) = median(Y(find(Y(:,dd)>Med),dd),'omitnan');
        IQR(dd) = 1.5 * (Q3(dd) - Q1(dd));
        L(dd) = Q1(dd) - IQR(dd);
        U(dd) = Q3(dd) + IQR(dd);
    end
    
    % Calculate 95% confidence interval
    CI95(dd) = confidenceInterval(Y(normdata{dd},dd),0.95,2); % 95%, 2 tailed
    % Calculate median
    MED(dd) = median(Y(normdata{dd},dd),'omitnan');
    % Calculate mean
    MEAN(dd) = mean(Y(normdata{dd},dd),'omitnan');
    
    % Calculate max/min
    try
        MAX(dd) = max(Y(normdata{dd},dd),[],'omitnan');
        MIN(dd) = min(Y(normdata{dd},dd),[],'omitnan');
    catch
        error(['Column ' num2str(dd) ' is not a distributive data set.']);
    end
    
    % Determine bin frequencies
    [HIST(:,dd), CENTER(:,dd)] = hist(Y(normdata{dd},dd),bins);
end

% Get maximum frequency for max width
mx = max(max(HIST));


%% CALCULATE PATCH COORDINATES
%=====================================

% FOR each column of data
for (dd = 1:size(Y,2))
    % Empty coordinate points
    xs = []; ys = [];
    % Max point
    xs = X(dd);
    ys = MAX(dd);
    
    % FOR each bin polygon (right side)
    for (bb = bins:-1:1)
        xs = [xs;    X(dd) + HIST(bb,dd)/mx*width/2];
        ys = [ys;    CENTER(bb,dd)];
    end
    
    % Min point
    xs = [xs; X(dd)];
    ys = [ys; MIN(dd)];
    
    % FOR each bin polygon (left side)
    for (bb = 1:bins)
        xs = [xs;    X(dd) - HIST(bb,dd)/mx*width/2];
        ys = [ys;    CENTER(bb,dd)];
    end
    
    % Connect the ends
    
    % Copy
    Xs(:,dd) = xs;
    Ys(:,dd) = ys;
end


%% DRAW PATCHES
%=====================================

H = ishold;
hold on;

% Violin Patch
%----------------------
% FOR each column of data
for (dd = 1:size(Y,2))
    % Draw patch
    patch(Xs(:,dd),Ys(:,dd),color{dd},'FaceAlpha',0.5);
        
    % 95% Condifence Interval Bars
    %----------------------
    UB(dd) = MEAN(dd) + CI95(dd);
    LB(dd) = MEAN(dd) - CI95(dd);
    w = 0.05 * width/2; % 5% max width
    patch([X(dd)-w X(dd)+w X(dd)+w X(dd)-w], [UB(dd) UB(dd) LB(dd) LB(dd)], 'k');
    pnt = X(dd) +w +w;
    
    
    % Mean Line
    %----------------------
    w = 0.125 * width/2; % 12.5% max width
    line([X(dd)-w*2 X(dd)+w*2], [MEAN(dd) MEAN(dd)],'LineWidth',2,'Color','k');

    % Median Point
    %----------------------
    plot(X(dd),MED(dd),'ow','MarkerSize',10);
    
    % Data Points
    %-----------------------
    if (IND_POINTS)
        plot(repmat(pnt,length(Y(normdata{dd},dd)),1), Y(normdata{dd},dd),'xk');
        if (~isempty(outliers{dd}))
            plot(repmat(pnt,length(Y(outliers{dd},dd)),1), Y(outliers{dd},dd),'xr');
        end
    end
    
    if (PLOT_BOX)
        % Draw Outlier Rectangle
        w = 0.11 * width/2; % 10% max width
        rectangle('Position',[X(dd)-w L(dd) w*2 (U(dd)-L(dd))],'EdgeColor','k');
    end
end

if (~H)
    hold off;
end


% Output statistics
if (nargout > 0)
    STATS.MEAN = MEAN;
    STATS.MED = MED;
    STATS.IQR = IQR;
    STATS.Q1 = Q1;
    STATS.Q3 = Q3;
    STATS.LQ = L;
    STATS.UQ = U;
    STATS.LB = LB;
    STATS.UB = UB;
end





end
