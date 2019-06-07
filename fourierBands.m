function [power, fbands] = fourierBands(sound, fs, varargin)
% function [power, fbands] = fourierBands(sound, fs, parameters...)
%
% Generate the magnitude spectrum for frequency bands and returns the
%      power in each band with the center frequencies.
%        Default: 1-octave widths
%
% Ex:  [power] = fourierBands(sound);
% Ex:  [power] = fourierBands(sound, 48000);
% Ex:  [power, fbands] = fourierBands(sound, 48000, 'bar','r');
% Ex:  [power, fbands] = fourierBands(sound, 48000, 'plot','-b');
% Ex:  [power, fbands] = fourierBands(sound, 48000, 'plot','-b','log');
%
% INPUTS:
% sound = Sound waveform (samples x 1)
% fs = sampling frequency (Default: 44100)
%
% Optional parameters:
%      'bar','S'   - Plots the spectral bands as a bar graph in current figure window,
%                       must be followed by color parameter, (ex: 'r').
%      'plot','S'  - Plots the spectral bands in current figure window,
%                       must be followed by plot parameters, (ex: 'b-o').
%      'log'       - Returns spectrum in Log magnitude (10*log10(x)) [dB FS].
%      'cuts',[L U] - Lower and upper frequency cut offs for band analysis.
%      'mean'      - Plots a mean power line to compare band outputs.
%
% Bandwidth parameter: [Choose one]
%      'octave'    - One-octave frequncy band width [Default].
%      'half'      - One-half octave frequncy bandwidth.
%      'third'     - One-third octave frequncy bandwidth.
%      'sixth'     - One-sixth octave frequncy bandwidth.
%      'band', B   - Custom bandwith in octaves, ex: '0.125' for 1/8 octave.
%
% Filter parameters:
%      'order', O  - Order of the filter [Default: 3]
%
% OUTPUTS:
% power = the relative power in each frequency band
% fbands = the frequency band center frequency
%
%   REQUIRES SIGNAL PROCESSING TOOLBOX
%
% @Author: Eugene Brandewie 08/04/2014
% Revised: 10/15/2016 - Added support for custom bandwidths and custom
%    filter orders. Removed unused 'norm' parameter.
% Revised: 01/11/2018 - Fixed error output for multi-column sound inputs.


%% Input Handling
%-----------------------------------------------------------
if (nargin < 1)
    error('Requires a waveform to analyze.');
elseif (nargin < 2)
    fs = 44100;
end

if (size(sound,2) > 1)
    error('Too many columns. Can only analyze one waveform channel.');
end



%% Handle optional input parameters
%-----------------------------------------------------------
LOG = 0;
PLOT = 0;
BAR = 0;
plotPars = '';
BAND = 1;
CUTS = [1 fs/2];
MEAN = 0;
ORDER = 3;

if nargin >= 3
    for ii = 1:nargin-2
        S = varargin{ii};
        if (~isnumeric(S))
            % LOG
            if strcmpi(S,'log')
                LOG = 1;
            end
            % PLOT
            if strcmpi(S,'plot')
                PLOT = 1;
                try
                    plotPars = varargin{ii+1};
                catch
                    error('Plotting parameters required...');
                end
            end
            % BAR
            if strcmpi(S,'bar')
                BAR = 1;
                try
                    plotPars = varargin{ii+1};
                catch
                    error('Bar plot parameters required...');
                end
            end
            % OCTAVE
            if strcmpi(S,'octave')
                BAND = 1;
            end
            % 1/2 OCTAVE
            if strcmpi(S,'half')
                BAND = 2;
            end
            % 1/3 OCTAVE
            if strcmpi(S,'third')
                BAND = 3;
            end
            % 1/6 OCTAVE
            if strcmpi(S,'sixth')
                BAND = 6;
            end
            % CUSTOM OCTAVES
            if strcmpi(S,'band')
                BAND = 1 / varargin{ii+1};
            end
            % FILTER ORDER
            if strcmpi(S,'order')
                ORDER = varargin{ii+1};
            end
            % CUTS
            if strcmpi(S,'cuts')
                CUTS = varargin{ii+1};
            end
            % MEAN
            if strcmpi(S,'mean')
                MEAN = 1;
            end
        end
    end
end



%% Band analysis
%-----------------------------------------------------------

if (BAND == 1) % Octave Band
    fCenter = (10.0^3) * ((2.0) .^ [-10:10]);
    fCenter = fCenter(fCenter>CUTS(1));
    fCenter = fCenter(fCenter<CUTS(2));
    fd = (2^(1/2));
    fUp = fCenter .* fd;
    fLow = fCenter ./ fd;
elseif (BAND == 2) % 1/2 Octave Band
    fCenter  = (10.0^3) * ((2.0) .^ ([-20:20]./2));
    fCenter = fCenter(fCenter>CUTS(1));
    fCenter = fCenter(fCenter<CUTS(2));
    fd = (2^(1/4));
    fUp = fCenter .* fd;
    fLow = fCenter ./ fd;
elseif (BAND == 3) % 1/3 Octave Band
    fCenter  = (10.0^3) * ((2.0) .^ ([-30:30]./3));
    fCenter = fCenter(fCenter>CUTS(1));
    fCenter = fCenter(fCenter<CUTS(2));
    fd = (2^(1/6));
    fUp = fCenter .* fd;
    fLow = fCenter ./ fd;
elseif (BAND == 6)% 1/6th Octave Band
    fCenter  = (10.0^3) * ((2.0) .^ ([-50:40]./6));
    fCenter = fCenter(fCenter>CUTS(1));
    fCenter = fCenter(fCenter<CUTS(2));
    fd = (2^(1/12));
    fUp = fCenter .* fd;
    fLow = fCenter ./ fd;
else % Custom Bandwith % TODO
    LP = BAND * log2(0.001); % 1 Hz
    HP = BAND * log2((fs/2)); % Nyquist    
    fCenter  = (10.0^3) * ((2.0) .^ ([LP:HP]./BAND));
    fCenter = fCenter(fCenter>CUTS(1));
    fCenter = fCenter(fCenter<CUTS(2));
    fd = (2^(1/(BAND*2)));
    fUp = fCenter .* fd;
    fLow = fCenter ./ fd;
end
nbands = length(fCenter);

% Nyquist frequency
nyq = fs/2;
% Prepare band power
bandPower = zeros(length(fCenter),1);
len = length(sound);

% FOR each center frequency
for ff = 1:length(fCenter)
    % Get the Nth-order butterworth filter
    [B,A] = butter(ORDER,[fLow(ff)/nyq,fUp(ff)/nyq]);
    % Filter waveform
    filtS = filter(B,A,sound);
    % Sum the power
    bandPower(ff,1) = sum(filtS.^2)/len;
end


%% IF logrithmic
%-----------------------------------------------------------
if (LOG)
    bandPower = 10*log10(bandPower);
    yLab = 'Power [dB FS]';
else
    yLab = 'Relative Power';
end


%% IF plotting
%-----------------------------------------------------------
if (PLOT || BAR)
    if (BAR)
        bar(bandPower,plotPars);
        set(gca,'XTick',[1:length(fCenter)]);
        set(gca,'XTickLabel',round(fCenter/5)*5);
        xLab = 'Frequency Band [Hz]';
    elseif (PLOT)
        plot(bandPower,plotPars);
        set(gca,'XTick',[1:length(fCenter)]);
        set(gca,'XTickLabel',round(fCenter/5)*5);
        xLab = 'Frequency [Hz]';
    end
    if (MEAN)
        line([0 length(fCenter)+1],[mean(bandPower) mean(bandPower)],'Color','k','LineWidth',2);
    end
    xlabel(xLab);
    ylabel(yLab);
end

% Output variables
power = bandPower;
fbands = fCenter;
























