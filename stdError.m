function [stdError] = stdError(matrix)
% function [stdError] = stdError(matrix)
% Calculates the standard error for each column in the given matrix.
%
% Eugene Brandewie 01/22/2013
% Updated: Sep 24, 2019 - Added 'omitnan' tag.

sd = std(matrix,'omitnan');
stdError = sd./sqrt(size(sd,1));