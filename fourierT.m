function [spectrum, freq, phase] = fourierT(sound, sF, P3, P4, P5, P6, P7)
% function [spectrum, freq, phase] = fourierT(sound, sF, parameters...)
%
% Generate the magnitude spectrum and phase angles of a waveform
%  using fast Fourier analysis.
%
% Ex:  [spect] = fourierT(SND)
% Ex:  [spect] = fourierT(SND, 32000)
% Ex:  [spect, f] = fourierT(SND, 'plot', '')
% Ex:  [spect, f] = fourierT(SND, 96000, 'plot', 'r', 1807)
%
% INPUTS:
% sound = Sound waveform (samples x 1)
% sF = sampling frequency (Default: 44100)
%
% Optional parameters:
%      'plot','S'  - Plots the magnitude spectrum in current figure window,
%                    must be followed by plot parameters, (ex: 'b-o')
%      'log'       - Returns spectrum in Log magnitude (20*Log10)
%      'N', S      - Performs an S-point fast fourier transform
%      'semi'      - Plots on semi-log X (frequency) axis
%
% OUTPUTS:
% spectrum = the fourier spectrum of the waveform
% freq = frequency points used for plotting the spectrum
% phase = phase value for each frequency point
%
%   REQUIRES SIGNAL PROCESSING TOOLBOX
%
% Eugene Brandewie 3/10/06
%  *Revised*   3/30/07 - Added phase ouput
%  *Revised*  11/19/07 - Added input parameters
%  *Revised*   2/22/08 - Added 'semi' input parameter



LOG = 0;
PLOT = 0;
SEMI = 0;
Np = 0;
N = length(sound);
if nargin >= 3
    for ii = 3:nargin
        eval(['S = P' num2str(ii) ';']);
        if (~isnumeric(S))
            if ~isempty(regexpi(S,'log'))
                LOG = 1;
            end
            if ~isempty(regexpi(S,'plot'))
                PLOT = 1;
                try
                    eval(['com = P' num2str(ii+1) ';']);
                catch
                    error('Plotting parameters required...');
                end
            end
            if ~isempty(regexpi(S,'N'))
                Np = 1;
                try
                    eval(['N = P' num2str(ii+1) ';']);
                catch
                    error('N-points required...');
                end
                
            end
            if ~isempty(regexpi(S,'semi'))
                SEMI = 1;
            end
        end
    end
end

if nargin < 2
    sF = 44100;
end

try
    SND = sound;
    dur = length(SND)/sF;
    npts=dur*sF;
    f=1:1/dur:sF/2;
    if N ~= length(sound)
        f=1:(sF)/N:sF/2;
    end
    warning off all;
    if Np
        X = fft(SND,N);
    else
        X = fft(SND);
    end
    spect=(abs(X));
    if LOG
        spect = 20*log10(abs(X));
    end
    phase = angle(X);
    warning on all;
catch
    error('ERROR: Cannot compute fourier transform...');
end

try
    sx = spect(2:round(npts/2));
    fx = f( 1:( length(f)-(length(f) - length(sx)) ) );
catch
    fx = f;
    sx = spect(2:length(f)+1);
end

try
    if PLOT
        if SEMI
            eval(['semilogx(fx, sx, ''' com ''');']);
        else
            eval(['plot(fx, sx, ''' com ''');']);
        end
        title('Frequency Content');
        xlabel('Frequency (Hz)');
        if LOG
            ylabel('Relative Magnitude (dB)');
        else
            ylabel('Relative Magnitude');
        end
    end
catch
    error('ERROR: Cannot plot the graph...');
end

spectrum = spect;
freq = f;
warning on;


