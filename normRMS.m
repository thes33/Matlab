function [newWave] = normRMS(waveform, rmsLevel, START, END);
%function [newWave] = normRMS(waveform, rmsLevel, START, END);
%
% Normalizes a waveform [columnar] to have the specified RMS amplitude.
% Can also specify the start and end point that is used for normalization 
% (such as the first 50 ms).  Multi-channel waveforms will have an 
% across-channel mean RMS equal to the given value, while maintaining 
% the ratio difference between channels.
%
% INPUTS:
%    waveform = Sound waveform (samples x channel) [columnar]
% 
% Optional parameters:
%    rmsLevel = Target across-channel mean RMS level. (Default: 0.5).
%    START  - Starting sample to take RMS calculation from.
%    END  - Ending sample to take RMS calculation from.
%  (Defaults to using all of the waveform)
%
% OUTPUTS:
%    newWave = the new scaled waveform
%
% Eugene Brandewie 01/22/2014

if (nargin < 2)
    rmsLevel = 0.5;
end

if (nargin < 3)
    START = 1;
    END = size(waveform,1);
end


% Get each channel RMS
curRMS = sqrt(mean(waveform(START:END,:).^2));
% Calculate ratio factor
Q = rmsLevel ./ mean(curRMS);
% Adjust amplitude based on mean RMS and scaling factor
newWave = waveform .* Q;













