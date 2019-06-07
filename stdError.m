function [stdError] = stdError(matrix)
% function [stdError] = stdError(matrix)
% Calculates the standard error for each column in the given matrix.
%
% Eugene Brandewie 01/22/2013

stdError = std(matrix)./sqrt(size(matrix,1));