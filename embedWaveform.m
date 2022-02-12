function [nX] = embedWaveform(X,Y,fs,varargin)
% function [nX] = embedWaveform(X,Y,fs,parameters...)
%  TODO: WIP!!!
%
% Embeds the waveform Y into the waveform X to create a new waveform, nX.
%
% INPUTS:
%   X = Sound waveform (samples x 1)
%   Y = Sound waveform to add (samples x 1)
%   fs = sampling frequency (Default: 48000)
%
% Optional parameters:
%      'SNR', N           - Apply the given SNR by increasing or reducing Y by the indicated decibel value in relation to X.
%      'onset', N (msec)  - Apply a starting offset in msec to Y in relation to the start of X.
%      'offset', N (msec) - Apply an ending offset in msec to X in relation to the end of Y.
%
% OUTPUTS:
%   nX = Modified X combined with Y
%
% @Author: Eugene Brandewie (Jan 4th, 2021)




%% Input Handling
%-----------------------------------------------------------
if (nargin < 2)
    error('Two waveforms required (X and Y).');
end

if (size(X,2) > size(X,1))
    error('X waveform must be columnar.');
end
if (size(Y,2) > size(Y,1))
    error('Y waveform must be columnar.');
end
if (size(X,2) > size(Y,2))
    error('X and Y must have same number of channels (columns).');
end

if (size(Y,1) > size(X,1))
    error('Y must have a shorter or equal duration to X.');
end


%% Handle optional input parameters
%-----------------------------------------------------------
SNR = 0;    % in dB
ONSET = 0; % In samples
OFFSET = 0; % In samples

if nargin >= 3
    for ii = 1:nargin-2
        S = varargin{ii};
        if (~isnumeric(S))
            % SNR
            if strcmpi(S,'SNR')
                SNR = varargin{ii+1};
            end
            % ONSET
            if strcmpi(S,'onset')
                ONSET = fs * varargin{ii+1}/1000;
                if (ONSET > length(X)), error('Onset cannot be greater than X waveform.'); end
            end
            % OFFSET
            if strcmpi(S,'offset')
                OFFSET = fs * varargin{ii+1}/1000;
            end
        end
    end
end

if (ONSET+size(Y,1)+OFFSET > length(X)), error('Additional Onset/Offset cannot be greater than X waveform.'); end



%% Apply Signal-to-Noise Ratio
%-----------------------------------------------------------
if (SNR ~= 0)    
    % Calculate the RMS value of X as baseline
    X_rms = mean(sqrt(mean(X.^2)));
    
    % Set Y RMS to match X RMS
    Y_rms = mean(sqrt(mean(Y.^2)));
    Qratio = X_rms ./ Y_rms;
    Yn = Y .* Qratio;
    
    % Apply Signal-to-Noise Ratio
    Y = Yn .* (10.^(SNR./20));    
end


%% Add Y to X
%-----------------------------------------------------------   













