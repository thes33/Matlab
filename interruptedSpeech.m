function [nWav, sqWave] = interruptedSpeech(wav, fs, interruptionRate, dutyCycle, startPhase)
% function [nWav] = interruptedSpeech(wav, fs, interruptionRate, dutyCycle)
%
%  Interrupted Speech Algorithm - Multiples a square wave with a waveform to create interrupted versions
%       of the signal. Returns the modified signal.
%
%   INPUTS:
%       wav - Waveform to modify.
%       fs - Sampling frequency of waveform.
%       interruptionRate - (Hz) [Default: 2) The frequency of the signal interruption.
%       dutyCycle - (Pecent) [Default: 50] The percent of the period in which the signal is present.
%       startPhase - [optional, Default: 0] The phase of the applied square wave.
%
%   OUTPUTS:
%       nWav - Modified waveform
%       sqWave - Square wave for analysis
%
%
% *REQUIRES SIGNAL PROCESSING TOOLBOX*
%

% Input Parameters
%======================================

if (nargin < 3)
    interruptionRate = 2.0;
end
if (nargin < 4)
    dutyCycle = 50;
end
if (nargin < 5)
    startPhase = 0;
end

% Square Wave Multiplication
%======================================
% Create time series
t = (1:length(wav))*1/fs; 
% Create square wave with input properties
sqWave= square((2*pi*interruptionRate*t) + deg2rad(startPhase), dutyCycle); 
% Remove negative values
sqWave(sqWave < 0) = 0; 
% Multiple with input waveform
nWav = sqWave .* wav';


