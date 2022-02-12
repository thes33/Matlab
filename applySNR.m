function [Xi] = applySNR(X,Y,SNR)
% function [Xi] = applySNR(X,Y,SNR)
%
% Applies the given RMS signal-to-noise ratio to X and Y by adjusting the RMS amplitude of X.
%
% INPUTS:
%   X = Sound waveform (samples x 1)
%   Y = Sound waveform (samples x 1)
%   SNR = Target SNR in decibels.
%
% OUTPUTS:
%   Xi = Modified X with adjusted amplitude.
%
% @Author: Eugene Brandewie Dec 18, 2020




%% Input Handling
%-----------------------------------------------------------
if (nargin < 3)
    error('Requires two waveforms to analyze and an SNR to apply.');
end

if (size(X,2) > 1) || (size(Y,2) > 1)
    error('Too many columns. Can only process one waveform channel.');
end



%% Apply Signal-to-Noise Ratio
%-----------------------------------------------------------

% Calculate the RMS value of Y as baseline
Y_rms = sqrt(mean(Y.^2));

% Set X RMS to match Y RMS
X_rms = mean(sqrt(mean(X.^2)));
Qratio = Y_rms ./ X_rms;
Xn = X .* Qratio;

% Apply Signal-to-Noise Ratio
Xi = Xn .* (10.^(SNR./20));
