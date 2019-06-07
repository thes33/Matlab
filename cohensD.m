function [D] = cohensD(X, varargin)
%function cohensD(X1, X2, X3 ...);
%
%  cohensD:  Calculates the Cohen's D effect size statistic between each given data column.
%
%     INPUTS:
%        X = Data to analyze, a column for each condition to compare. Must be at least two columns.
%           Data can be a single RxC matrix, or each data column can be inputted separately (for varying N).
%
%     OUTPUTS:
%        D = CxC N - Cohen's D statistic for each column combination.
%
%  @Author: Eugene Brandewie - April 25, 2019
%
%  VERSION:
%   1.0 - Initial Release
%



%% PROCESS INPUTS
%===================================
% Set defaults
CELL = false;

% Ensure enough columns
if (nargin < 2 && size(X,2) < 2)
    error('Not enough input arguments: need at least 2 data columns.');
end

% IF single matrix
%----------------------------
if (size(X,2) > 1)
    
    % Create cell array
    DATA = {};
    % FOR each additional column
    for (ii = 1:size(X,2))
        DATA(ii) = {X(:,ii)};
    end
    
else % ELSE separate columns
    %----------------------------
    % Create cell array
    DATA = {};
    DATA(1) = {X};
    
    % FOR each additional column
    if nargin > 1
        for ii = 1:(nargin-1)
            Q = varargin{ii};
            % Data columns
            if (isnumeric(Q))
                DATA(ii+1) = {Q};
            end
        end
    end
end



%% CALCULATE STATISTICS
%=====================================

% FOR each column of data
for (dd = 1:length(DATA))
    col1 = dd;
    
    % FOR each column of data
    for (cc = 1:length(DATA))
        col2 = cc;
        
        % Ignore same column
        if (col1 == col2)
            D(col1,col2) = 0; continue;
        end
        clear M1 M2 S1 S2 N1 N2 Sp;
                
        % Means
        M1 = mean(DATA{col1},'omitnan');
        M2 = mean(DATA{col2},'omitnan');
        
        % Standard Deviations
        S1 = std(DATA{col1},'omitnan');
        S2 = std(DATA{col2},'omitnan');
                
        % Number
        N1 = length(DATA{col1}) - sum(isnan(DATA{col1}));
        N2 = length(DATA{col2}) - sum(isnan(DATA{col2}));
        
        % Pooled Standard Deviation
        Sp = sqrt( (((N1-1)*S1*S1) + ((N2-1)*S2*S2)) / (N1+N2));
             
        % Cohen's D
        D(col1,col2) = (M1 - M2) / Sp;
    end
end

































