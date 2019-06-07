function [RMS] = calcRMS(signal)
% function [RMS] = calcRMS(signal)
% Calculates the RMS level of the signal.
%
% Eugene Brandewie 10/1/07

N = signal.^2;
N = mean(N);
N = sqrt(N);

RMS = N;
