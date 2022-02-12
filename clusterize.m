function [Xh] = clusterize(X, Xvar, type)
% function [Xh] = clusterize(X, Xvar, type)
%
% Add variance around values to create a cluster of data points
%   and increase readability when graphing.
%
% INPUTS:
%    X = X values to add variance.
%    Xvar = Max amount of +/- adjustments applied to X values.
%
% Optional parameters:
%      type = 'random' or 'uniform' - Choose between a 'random' variance or a uniform [default] cluster.
%
% OUTPUTS:
%    Xh = Modified values of X
%
% @Author: Eugene Brandewie - Dec 7th, 2021



%% INPUT HANDLING
%=========================================

if (nargin < 3)
    type = 'uniform';
end


%% Apply random variance
%=========================================

% IF, uniform variance
if (strcmp(type,'uniform'))
    Xe = linspace(-Xvar,Xvar,length(X));
    for (xx = 1:length(X))
        Xh(xx) = X(xx) + Xe(xx);
    end
    
else % ELSE, random variance
    for (xx = 1:length(X))
        v = (rand()*Xvar*2) - Xvar;
        Xh(xx) = X(xx) + v;
    end
end

