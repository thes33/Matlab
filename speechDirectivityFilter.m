function [modSpeechWav] = speechDirectivityFilter(speechWav, fs, angle);
% function [modSpeechWav] = speechDirectivityFilter(speechWav);
%
%  Applies a directivty filter to the given waveform. This filter
%   is based on the normal speech directivity filter presented
%   in Monson et al (2012).
%
% INPUT:
%    speechWav - Waveform of the speech stimulus to process, 
%                 single-channel columnar.
%    fs - Sampling rate of waveform
%    angle - Azimuth angle to simulate with the directivity filter.
%
% References:
% -------------
% Monson, B. B., Hunter, E. J., & Story, B. H. (2012). 
%   Horizontal directivity of low- and high-frequency 
%   energy in speech and singing. JASA, 132(1), 433–441.

% Set up table for normal speech adjustments
levelAdjust = [0,0,0,0,0,0,0,0;0.800000000000000,0.700000000000000,0.600000000000000,-0.200000000000000,-1, ...
    -1.70000000000000,-1,-0.100000000000000;-0.100000000000000,-0.200000000000000,0.100000000000000, ...
    -0.700000000000000,-0.100000000000000,-0.700000000000000,-1.70000000000000,-0.500000000000000;-0.400000000000000, ...
    -0.500000000000000,0,-1.40000000000000,-0.800000000000000,-1.60000000000000,-3.10000000000000,-2.30000000000000; ...
    -0.600000000000000,-0.800000000000000,-0.100000000000000,-1.50000000000000,-2.50000000000000,-3.90000000000000, ...
    -5.80000000000000,-4.40000000000000;-1.10000000000000,-1.60000000000000,-0.700000000000000,-1.40000000000000, ...
    -4.70000000000000,-5.30000000000000,-7,-5.10000000000000;-1.60000000000000,-2.30000000000000,-1.60000000000000, ...
    -1.60000000000000,-7,-6.30000000000000,-8.70000000000000,-7.20000000000000;-2.10000000000000,-3.10000000000000, ...
    -2.80000000000000,-2.50000000000000,-7.60000000000000,-9,-11.5000000000000,-10.1000000000000;-2.70000000000000, ...
    -3.90000000000000,-4.40000000000000,-4.60000000000000,-7.90000000000000,-12.7000000000000,-15,-14.2000000000000; ...
    -3.10000000000000,-4.50000000000000,-5.50000000000000,-7,-9.30000000000000,-14.5000000000000,-17.8000000000000, ...
    -17.8000000000000;-3.50000000000000,-4.90000000000000,-6,-7.90000000000000,-13,-16.7000000000000,-21.5000000000000, ...
    -21.7000000000000;-3.60000000000000,-5,-6,-6.90000000000000,-15,-21,-25,-23.5000000000000;-3.60000000000000, ...
    -5.10000000000000,-5.90000000000000,-6.30000000000000,-13.1000000000000,-18.9000000000000,-26.3000000000000, ...
    -26.7000000000000];
angles = [0;15;30;45;60;75;90;105;120;135;150;165;180];
freqs = [125,250,500,1000,2000,4000,8000,16000];

% Use interpolation to find values for given angle
for (ff = 1:length(freqs))
    levelAdj(ff) = interp1(angles,levelAdjust(:,ff),angle,'pchip');
end

% Create filter centers
fd = (2^(1/2));
fUp = freqs .* fd;
fLow = freqs ./ fd;
nyq = fs/2;

% FOR each center frequency
for ff = 1:length(freqs)
    % Get the Nth-order butterworth filter
    [B,A] = butter(3,[fLow(ff)/nyq,fUp(ff)/nyq]);
    % Filter waveform using octave filter
    filteredWav(:,ff) = filter(B,A,speechWav);
    % Apply associated power changes to each filter channel
    filteredWav(:,ff) = filteredWav(:,ff) .* (10.^(levelAdj(ff)/20));    
end

% Return the filtered speech waveform
modSpeechWav = sum(filteredWav')';














