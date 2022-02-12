% Copyright (C) 2005, 2006 William Poetra Yoga Hadisoeseno
% Copyright (C) 2011 Nir Krakauer
%
% This program is free software; you can redistribute it and/or modify it under
% the terms of the GNU General Public License as published by the Free Software
% Foundation; either version 3 of the License, or (at your option) any later
% version.
%
% This program is distributed in the hope that it will be useful, but WITHOUT
% ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
% FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
% details.
%
% You should have received a copy of the GNU General Public License along with
% this program; if not, see <http://www.gnu.org/licenses/>.

% -*- texinfo -*-
% @deftypefn {Function File} {[@var{b}, @var{bint}, @var{r}, @var{rint}, @var{stats}] =} regress (@var{y}, @var{X}, [@var{alpha}])
% Multiple Linear Regression using Least Squares Fit of @var{y} on @var{X}
% with the model @code{y = X * beta + e}.
%
% Here,
%
% @itemize
% @item
% @code{y} is a column vector of observed values
% @item
% @code{X} is a matrix of regressors, with the first column filled with
% the constant value 1
% @item
% @code{beta} is a column vector of regression parameters
% @item
% @code{e} is a column vector of random errors
% @end itemize
%
% Arguments are
%
% @itemize
% @item
% @var{y} is the @code{y} in the model
% @item
% @var{X} is the @code{X} in the model
% @item
% @var{alpha} is the significance level used to calculate the confidence
% intervals @var{bint} and @var{rint} (see `Return values' below). If not
% specified, ALPHA defaults to 0.05
% @end itemize
%
% Return values are
%
% @itemize
% @item
% @var{b} is the @code{beta} in the model
% @item
% @var{bint} is the confidence interval for @var{b}
% @item
% @var{r} is a column vector of residuals
% @item
% @var{rint} is the confidence interval for @var{r}
% @item
% @var{stats} is a row vector containing:
%
%   @itemize
%   @item The R^2 statistic
%   @item The F statistic
%   @item The p value for the full model
%   @item The estimated error variance
%   @end itemize
% @end itemize
%
% @var{r} and @var{rint} can be passed to @code{rcoplot} to visualize
% the residual intervals and identify outliers.
%
% NaN values in @var{y} and @var{X} are removed before calculation begins.
%
% @end deftypefn

% References:
% - Matlab 7.0 documentation (pdf)
% - 《大学数学实验》 姜启源 等 (textbook)
% - http://www.netnam.vn/unescocourse/statistics/12_5.htm
% - wsolve.m in octave-forge
% - http://www.stanford.edu/class/ee263/ls_ln_matlab.pdf

function [b, bint, r, rint, stats] = regress (y, X, alpha)

if (nargin < 2 || nargin > 3)
    print_usage;
end

if (~ ismatrix (y))
    error ("regress: y must be a numeric matrix");
end
if (~ ismatrix (X))
    error ("regress: X must be a numeric matrix");
end

if (size(y,2) ~= 1)
    error ("regress: y must be a column vector");
end

if (size(y,1) ~= size(X,1))
    error ("regress: y and X must contain the same number of rows");
end

if (nargin < 3)
    alpha = 0.05;
elseif (~ isscalar (alpha))
    error ("regress: alpha must be a scalar value")
end

notnans = ~ logical (sum (isnan ([y X]), 2));
y = y(notnans);
X = X(notnans,:);

[Xq Xr] = qr (X, 0);
pinv_X = Xr \ Xq';

b = pinv_X * y;

if (nargout > 1)
    
    n = rows (X);
    p = columns (X);
    dof = n - p;
    t_alpha_2 = tinv (alpha / 2, dof);
    
    r = y - X * b; % added -- Nir
    SSE = sum (r .^ 2);
    v = SSE / dof;
    
    % c = diag(inv (X' * X)) using (economy) QR decomposition
    % which means that we only have to use Xr
    c = diag (inv (Xr' * Xr));
    
    db = t_alpha_2 * sqrt (v * c);
    
    bint = [b + db, b - db];
    
end

if (nargout > 3)
    
    dof1 = n - p - 1;
    h = sum(X.*pinv_X', 2); %added -- Nir (same as diag(X*pinv_X), without doing the matrix multiply)
    
    % From Matlab's documentation on Multiple Linear Regression,
    %   sigmaihat2 = norm (r) ^ 2 / dof1 - r .^ 2 / (dof1 * (1 - h));
    %   dr = -tinv (1 - alpha / 2, dof) * sqrt (sigmaihat2 .* (1 - h));
    % Substitute
    %   norm (r) ^ 2 == sum (r .^ 2) == SSE
    %   -tinv (1 - alpha / 2, dof) == tinv (alpha / 2, dof) == t_alpha_2
    % We get
    %   sigmaihat2 = (SSE - r .^ 2 / (1 - h)) / dof1;
    %   dr = t_alpha_2 * sqrt (sigmaihat2 .* (1 - h));
    % Combine, we get
    %   dr = t_alpha_2 * sqrt ((SSE * (1 - h) - (r .^ 2)) / dof1);
    
    dr = t_alpha_2 * sqrt ((SSE * (1 - h) - (r .^ 2)) / dof1);
    
    rint = [r + dr, r - dr];
    
end

if (nargout > 4)
    
    R2 = 1 - SSE / sum ((y - mean (y)) .^ 2);
    %    F = (R2 / (p - 1)) / ((1 - R2) / dof);
    F = dof / (p - 1) / (1 / R2 - 1);
    pval = 1 - fcdf (F, p - 1, dof);
    
    stats = [R2 F pval v];
    
end

end
