function [err, lb, ub, df] = confidenceInterval(matrix, interval, tails)
% function [err, lb, ub, df] = confidenceInterval(matrix, interval, tails)
%
% Calculates the confidence interval for the given matrix. Uses Student's
%   t distribution for small samples (< 30), uses Z-score distribution for
%   larger samples. Processes each column of the input matrix.
%
% INPUT:
%  matrix - Matrix of data points to analyze.
%  interval - (optional) Confidence interval [Default: 0.95, 95%].
%     Options: One-tailed: 0.85, 0.90, 0.95, 0.975, 0.99, 0.995, 0.9975,
%                          0.999, 0.9995.
%              Two-tailed: 0.70, 0.80, 0.90, 0.95, 0.98, 0.99, 0.995,
%                          0.9975, 0.999.
%  tails - (optional) 1 or 2 tailed [Default: 2].
%
% OUTPUT:
%  err - Errors (in 'matrix' units).
%  lb - Lower bounds (mean - err).
%  ub - Upper bounds (mean + err).
%  df - Degrees of freedom.
%
% @Author: Eugene Brandewie 04/05/2013
% @Email: eugene.brandewie@gmail.com
% @Version: 1.3.1



% VERSION HISTORY
%------------------------------------------
% 1.1: 07/20/2016 - No longer requires Statistics Toolbox. Uses a
%    look-up table for T and Z distributions.
% 1.2: 10/03/2016 - Now ignores NaN values like nanmean.
% 1.3: 10/07/2016 - Processes input matrix columns separately.
% 1.3.1: 10/10/2016 - Fixed 1-D row-arrays to transpose before processing.




%% SET DEFAULT PARAMETERS
%-----------------------------------------
if (nargin < 1)
    error('Must supply a matrix to process.');
end
if (nargin < 2)
    interval = 0.95;
end
if (nargin < 3)
    tails = 2;
end


% Check NaN values
if (isempty(matrix(~isnan(matrix))))
    warning(['All matrix values are NaN.' ...
            ' No calculable confidence interval.']);
    err = NaN;
    lb = NaN;
    ub = NaN;
    df = 0;
    return;
end


% Confidence Interval options
CI1 = [0.85, 0.90, 0.95, 0.975, 0.99, 0.995, 0.9975, 0.9990, 0.9995];
CI2 = [0.70, 0.80, 0.90, 0.950, 0.98, 0.990, 0.9950, 0.9975, 0.9990];

% Check interval within parameters
if (tails==1)
    if (sum(CI1==interval) == 0)
        error(['''interval'' must be one of the following for one-tailed tests:  ' ...
            num2str(unique(CI1))]);
    end
elseif (tails==2)
    if (sum(CI2==interval) == 0)
        error(['''interval'' must be one of the following for two-tailed tests:  ' ...
            num2str(unique(CI2))]);
    end
end


% Transpose if 1-D row-array
if (size(matrix,1) == 1)
    matrix = matrix';
end


%% PROCESS EACH COLUMN
%-----------------------------------------

% FOR each matrix column
for (cc = 1:size(matrix,2))
    % Get working column
    workCol = matrix(:,cc);
    
    % Remove NaN values
    workCol = workCol(~isnan(workCol));
    if (isempty(workCol))
        warning(['All values are NaN for column ' num2str(cc) '.']);
    end
    
    % CHECK SAMPLE SIZE AND INTERVAL
    %-----------------------------------------

    % Get sample size
    sz = size(workCol,1);
    
    % Get degrees of freedom
    df(cc) = sz-1;
    if (df(cc) == 0)
        warning(['Column ' num2str(cc) ' has only one data point.' ...
            ' No calculable confidence interval.']);
        err(cc) = 0.0;
        lb(cc) = mean(workCol);
        ub(cc) = mean(workCol);
        df(cc) = 0;
        continue;
    end
    
    % Check DF
    if (df(cc) > 30)
        USE_Z = true;
    else
        USE_Z = false;
    end
    
    
    % Find index of confidence interval
    if (tails==1)
        ci = find(CI1==interval);
        ci = ci(1);
    elseif (tails==2)
        ci = find(CI2==interval);
        ci = ci(1);
    end
    
    
    % BASIC STATISTICS
    %-----------------------------------------
    
    % Get Mean
    mn = mean(workCol);
    % Get Standard deviation
    sd = std(workCol);
    % Get Standard error
    se = sd./sqrt(sz);
    
    
    % USE T DISTRIBUTION
    %--------------------------------------------------
    if (USE_Z == false)
        % Degrees of Freedom
        DF = [1;2;3;4;5;6;7;8;9;10;11;12;13;14;15;16;17;18;19;20;21;22;23;24;25; ...
            26;27;28;29;30];
        % Student's T Look-up Table
        TDIST = [1.963,3.078,6.314,12.710,31.820,63.660,127.300,318.300,636.600; ...
            1.386,1.886,2.920,4.303,6.965,9.925,14.090,22.330,31.600;
            1.250,1.638,2.353,3.182,4.541,5.841,7.453,10.210,12.920;
            1.190,1.533,2.132,2.776,3.747,4.604,5.598,7.173,8.610;
            1.156,1.476,2.015,2.571,3.365,4.032,4.773,5.893,6.869;
            1.134,1.440,1.943,2.447,3.143,3.707,4.317,5.208,5.959;
            1.119,1.415,1.895,2.365,2.998,3.499,4.029,4.785,5.408;
            1.108,1.397,1.860,2.306,2.896,3.355,3.833,4.501,5.041;
            1.100,1.383,1.833,2.262,2.821,3.250,3.690,4.297,4.781;
            1.093,1.372,1.812,2.228,2.764,3.169,3.581,4.144,4.587;
            1.088,1.363,1.796,2.201,2.718,3.106,3.497,4.025,4.437;
            1.083,1.356,1.782,2.179,2.681,3.055,3.428,3.930,4.318;
            1.079,1.350,1.771,2.160,2.650,3.012,3.372,3.852,4.221;
            1.076,1.345,1.761,2.145,2.624,2.977,3.326,3.787,4.140;
            1.074,1.341,1.753,2.131,2.602,2.947,3.286,3.733,4.073;
            1.071,1.337,1.746,2.120,2.583,2.921,3.252,3.686,4.015;
            1.069,1.333,1.740,2.110,2.567,2.898,3.222,3.646,3.965;
            1.067,1.330,1.734,2.101,2.552,2.878,3.197,3.610,3.922;
            1.066,1.328,1.729,2.093,2.539,2.861,3.174,3.579,3.883;
            1.064,1.325,1.725,2.086,2.528,2.845,3.153,3.552,3.850;
            1.063,1.323,1.721,2.080,2.518,2.831,3.135,3.527,3.819;
            1.061,1.321,1.717,2.074,2.508,2.819,3.119,3.505,3.792;
            1.060,1.319,1.714,2.069,2.500,2.807,3.104,3.485,3.767;
            1.059,1.318,1.711,2.064,2.492,2.797,3.091,3.467,3.745;
            1.058,1.316,1.708,2.060,2.485,2.787,3.078,3.450,3.725;
            1.058,1.315,1.706,2.056,2.479,2.779,3.067,3.435,3.707;
            1.057,1.314,1.703,2.052,2.473,2.771,3.057,3.421,3.690;
            1.056,1.313,1.701,2.048,2.467,2.763,3.047,3.408,3.674;
            1.055,1.311,1.699,2.045,2.462,2.756,3.038,3.396,3.659;
            1.055,1.310,1.697,2.042,2.457,2.750,3.030,3.385,3.646;];
        
        % Look-up t value
        t = TDIST(df(cc),ci);
        % Get error adjustment
        err(cc) = se * t;
        
    else
        % ELSE USE Z DISTRIBUTION
        %--------------------------------------------------
        % Z-Score Look-up Table
        ZDIST = [1.036,1.282,1.645,1.960,2.326,2.576,2.807,3.090,3.291];
        % Look-up Z value
        z = ZDIST(ci);
        % Get error adjustment
        err(cc) = se * z;
        
    end
    
    
    % CALCULATE BOUNDS
    %-----------------------------------------
    % Compute lower-bound
    lb(cc) = mn - err(cc);
    % Compute upper-bound
    ub(cc) = mn + err(cc);
    
end % END for each column















