function [STATS] = boxPlot(X, Y, varargin);
%function boxPlot(X, Y, parameters...);
%
%  boxPlot:  Plot box plots showing quartiles and ranges of data sets.
%
%     INPUTS:
%        X = x coordinates for each column (center of violins).
%        Y = Y data value columns to determine means and errors.
%
%     OPTIONAL INPUTS:
%        'width', N = Desired maximum width of the violin plot densities. [Default 0.8].
%        'points' = Plots individual data points as X's. Outliers in red.
%        'labels', C = Add labels to list of individual data points. Cell array of strings.
%        'nobox' = Does not plot interquartile range box.
%        'noci' = Does not plot 95% confidence intervals.
%        'iqr' = Plots the extra-interquartile range box (1.5 x IQR) for outlier detection.
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
%  @Author: Eugene Brandewie Sep 28, 2018
%
%  VERSION:
%   1.0 - Initial Release
%   1.1 - Oct 3,2018 - Fixed issues with uniform and near-uniform data sets.
%   1.2 - Mar 8,2019 - Added support for data point labels.
%   1.2.1 - Mar 20,2019 - Fixed issue with labels and multiple columns.
%   1.2.2 - May 9,2019 - Fixed issue with points and labels.
%


%% PROCESS INPUTS
%===================================
% Set defaults
color = 'k';
width = 0.8;
bins = 5;
IND_POINTS = false;
DATA_LABELS = {};
IGNORE_OUTLIERS = false;
PLOT_IQR = false;
BOX = true;
SHOW_CI = true;

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
            % Width
            if (strcmp(Q,'width'))
                width = varargin{ii+1};
            end
            % No IQR Box
            if (strcmp(Q,'nobox'))
                BOX = false;
            end
            % No CI95%
            if (strcmp(Q,'noci'))
                SHOW_CI = false;
            end
            % Individual Data Points
            if (strcmp(Q,'points'))
                IND_POINTS = true;
            end
            % Individual Data Points
            if (strcmp(Q,'labels'))
                DATA_LABELS = varargin{ii+1};
                if (~iscell(DATA_LABELS))
                    error('''labels'' must be a cell array of strings.');
                end
                if (length(DATA_LABELS) ~= size(Y,1))
                    error('Data labels must be same length as number of points.');
                end
                IND_POINTS = true;
            end
            % Ignore Outliers
            if (strcmp(Q,'outliers'))
                IGNORE_OUTLIERS = true;
            end
            % IQR
            if (strcmp(Q,'iqr'))
                PLOT_IQR = true;
            end
        end
    end
end




%% CALCULATE STATISTICS
%=====================================

% FOR each column of data
for (dd = 1:size(Y,2))
    
    % PROCESS & REMOVE  OUTLIERS
    %-------------------------
    if (~IGNORE_OUTLIERS)
        % Determine Outliers
        Med = median(Y(:,dd),'omitnan');
        Q1(dd) = median(Y(find(Y(:,dd)<Med),dd),'omitnan');
        Q3(dd) = median(Y(find(Y(:,dd)>Med),dd),'omitnan');
        
        % Deal with non distributive sets
        if (isnan(Q1(dd)) || isnan(Q3(dd)))
            ySort = sort(Y(:,dd));
            % IF even
            if  (mod(length(ySort),2) == 0)
                Q1(dd) = median(ySort(1:length(ySort)/2),'omitnan');
                Q3(dd) = median(ySort((length(ySort)/2)+1:end),'omitnan');
            else
                Q1(dd) = median(ySort(1:(length(ySort)/2)-1),'omitnan');
                Q3(dd) = median(ySort((length(ySort)/2)+1:end),'omitnan');
            end
        end
        
        IQR(dd) = 1.5 * (Q3(dd) - Q1(dd));
        L(dd) = Q1(dd) - IQR(dd);
        U(dd) = Q3(dd) + IQR(dd);
        
        outIndx = find(Y(:,dd) > U(dd) | Y(:,dd) < L(dd));
        outliers(dd) = {outIndx};
        if (~isempty(outliers{dd}))
            warning([num2str(length(outliers{dd})) ' outliers found in column ' num2str(dd) '.']);
        end
        normIndx = find(Y(:,dd) < U(dd) & Y(:,dd) > L(dd));
        normdata(dd) = {normIndx};
        
        % Data point labels
        if (~isempty(DATA_LABELS))
            NORM_LABELS(dd) = {DATA_LABELS(normIndx)};
            OUT_LABELS(dd) = {DATA_LABELS(outIndx)};
        end
        
        % If all data outlier, distribution is low, so include all
        if (isempty(normdata{dd}))
            normdata(dd) = {1:length(Y(:,dd))};
            warning(['Column ' num2str(dd) ' is not a distributive data set.']);
            NORM_LABELS(dd) = {DATA_LABELS(1:length(Y(:,dd)))};
        end
        
        % IGNORING OUTLIERS
        %-------------------------
    else % Else ignore outliers
        Med = median(Y(:,dd),'omitnan');
        normdata(dd) = {1:length(Y(:,dd))};
        outliers(dd) = {[]};
        
        % Data point labels
        if (~isempty(DATA_LABELS))
            NORM_LABELS(dd) = {DATA_LABELS(1:length(Y(:,dd)))};
        end
        
        Q1(dd) = median(Y(find(Y(:,dd)<Med),dd),'omitnan');
        Q3(dd) = median(Y(find(Y(:,dd)>Med),dd),'omitnan');
        
        % Deal with non distributive sets
        if (isnan(Q1(dd)) || isnan(Q3(dd)))
            ySort = sort(Y(:,dd));
            % IF even
            if  (mod(length(ySort),2) == 0)
                Q1(dd) = median(ySort(1:length(ySort)/2),'omitnan');
                Q3(dd) = median(ySort((length(ySort)/2)+1:end),'omitnan');
            else
                Q1(dd) = median(ySort(1:(length(ySort)/2)-1),'omitnan');
                Q3(dd) = median(ySort((length(ySort)/2)+1:end),'omitnan');
            end
        end
        
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
    MAX(dd) = max(Y(normdata{dd},dd),[],'omitnan');
    MIN(dd) = min(Y(normdata{dd},dd),[],'omitnan');
end



%% DRAW PATCHES
%=====================================

H = ishold;
hold on;

% FOR each column of data
for (dd = 1:size(Y,2))
    
    % Center Line
    %----------------------
    line([X(dd) X(dd)], [MIN(dd) MAX(dd)], 'Color','k');
    
    % Min/Max Lines
    %----------------------
    w = 0.65 * width/2; % 5% max width
    line([X(dd)-w X(dd)+w], [MIN(dd) MIN(dd)], 'Color','k');
    line([X(dd)-w X(dd)+w], [MAX(dd) MAX(dd)], 'Color','k');
    
    % Interquartile Box
    %----------------------
    if (BOX)
        w = width/2;
        rectangle('Position',[X(dd)-w Q1(dd) width Q3(dd)-Q1(dd)],'EdgeColor','k','LineWidth',1);
    end
    
    
    % 1.5x Interquartile Box
    %----------------------
    if (PLOT_IQR)
        w = 0.65 * width/2;
        rectangle('Position',[X(dd)-w L(dd) width*0.65 U(dd)-L(dd)],'Curvature',1,'EdgeColor',[0.7 0.7 0.7],'LineWidth',1);
    end
    
    % Median Line
    %----------------------
    w = width/2;
    line([X(dd)-w X(dd)+w], [MED(dd) MED(dd)], 'Color','k');
    
    
    % 95% Condifence Interval Bars
    %----------------------
    if (SHOW_CI)
        UB(dd) = MEAN(dd) + CI95(dd);
        LB(dd) = MEAN(dd) - CI95(dd);
        w = 0.06 * width/2; % 5% max width
        patch([X(dd)-w X(dd)+w X(dd)+w X(dd)-w], [UB(dd) UB(dd) LB(dd) LB(dd)], 'k');
    end
    
    % Mean Line
    %----------------------
    w = 0.125 * width/2; % 12.5% max width
    line([X(dd)-w*2 X(dd)+w*2], [MEAN(dd) MEAN(dd)],'LineWidth',2,'Color','k');
    
    % Data Points
    %-----------------------
    if (IND_POINTS)
        pnt = X(dd) +w +w;
        
        if (exist('NORM_LABELS','var') && ~isempty(NORM_LABELS{dd}))
            text(repmat(pnt,length(Y(normdata{dd},dd)),1),Y(normdata{dd},dd),NORM_LABELS{dd});
        else
            plot(repmat(pnt,length(Y(normdata{dd},dd)),1), Y(normdata{dd},dd),'xk');
        end
        
        if (~isempty(outliers{dd}))
            if (exist('OUT_LABELS','var') && ~isempty(OUT_LABELS{dd}))
                text(repmat(pnt,length(Y(outliers{dd},dd)),1),Y(outliers{dd},dd),OUT_LABELS{dd},'Color','r');
            else
                plot(repmat(pnt,length(Y(outliers{dd},dd)),1), Y(outliers{dd},dd),'xr');
            end
        end
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
