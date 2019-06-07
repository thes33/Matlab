function [fsound, freqs, summed] = filterBank(sound, fs, varargin)
% function [fsound, freqs, summed] = filterBank(sound, fs, parameters...)
%
% Filters the given sound waveform through a bank of filters
%      and returns an array of the filtered waveforms.
%
% INPUTS:
% sound = Sound waveform (samples x 1)
% fs = sampling frequency (Default: 44100)
%
% Bandwidth parameter: [Choose one]
%      'octave'    - One-octave frequncy band width [Default].
%      'half'      - One-half octave frequncy bandwidth.
%      'third'     - One-third octave frequncy bandwidth.
%      'sixth'     - One-sixth octave frequncy bandwidth.
%      'band', B   - Custom bandwith in octaves, ex: '0.125' for 1/8 octave.
%      'erb', C    - Bandwidths based on ERBs (equivalent rectangular
%                      bandwidths). Using Moore and Glasberg (1990) formula.
%                      C - Number of ERB channels. Use CUTS to narrow range.
%
% Optional filter parameters:
%      'order', O   - Order of the Butterworth filters [Default: 3]
%      'cuts',[L U] - Lower and upper frequency cut-offs for band analysis.
%                       [Default: 50 Hz to 20 kHz]
%      'type', T    - Choose the filter type: 'butter'  for Butterworth [Default]
%                       or 'gamma' for a gammatone filterbank.
%      'align'      - Phase aligns the fine-structure and envelope for 
%                       gamma-tone filtering such that peak occurs at t=0.
%                       [Default: off]
%
% OUTPUTS:
% fsound(:,N) = N filtered waveforms, one for each band.
% freqs = Center frequencies of each filter's band.
% summed = (optional) Summed waveform of all channels combined.
%
%   REQUIRES SIGNAL PROCESSING TOOLBOX
%
% @Author: Eugene Brandewie 12/15/2016
% Modified: EJB 3/10/2017: Added ERB bandwidths and filter types.




%% REQUIRED INPUTS
%====================================================
if (nargin < 1)
    error('Requires a waveform to filter.');
elseif (nargin < 2)
    fs = 44100;
end



%% OPTIONAL INPUT PARAMETERS
%====================================================
BAND = 1;               % Ratio of octave for bandwidths
CUTS = [50 20000];      % Frequency limits
ORDER = 3;              % Filter-order
ERB = 0;                % ERB flag
CHANNELS = 1;           % Number of ERB channels
FILTTYPE = 'butter';    % Filter Type: 'butter' or 'gamma'
PHASEALIGN = false;     % Phase-alignment for gamma-tone

if nargin > 2
    for ii = 1:nargin-2
        S = varargin{ii};
        if (~isnumeric(S))
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
            % ERB
            if strcmpi(S,'erb')
                BAND = NaN;
                ERB = 1;
                CHANNELS = varargin{ii+1};
            end
            % FILTER ORDER
            if strcmpi(S,'order')
                ORDER = varargin{ii+1};
            end
            % CUT
            if strcmpi(S,'cuts')
                CUTS = varargin{ii+1};
            end
            % ALIGN
            if strcmpi(S,'align')
                PHASEALIGN = true;
            end
        end
    end
end



%% BAND CREATION
%====================================================

if (BAND == 1) % Octave Band
    fCenter = (10.0^3) * ((2.0) .^ [-10:10]);
    fCenter = fCenter(fCenter>CUTS(1));
    fCenter = fCenter(fCenter<CUTS(2));
    fd = (2^(1/2));
    fUp = fCenter .* fd;
    fLow = fCenter ./ fd;
    
elseif (BAND == 3) % 1/2 Octave Band
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
    
elseif (~ERB) % Custom Octave-based Bandwith
    LP = BAND * log2(0.001); % 1 Hz
    HP = BAND * log2((fs/2)); % Nyquist
    fCenter  = (10.0^3) * ((2.0) .^ ([LP:HP]./BAND));
    fCenter = fCenter(fCenter>CUTS(1));
    fCenter = fCenter(fCenter<CUTS(2));
    fd = (2^(1/(BAND*2)));
    fUp = fCenter .* fd;
    fLow = fCenter ./ fd;
    
elseif (ERB) % ERB Bandwidths
    % Calculate max ERB filter number
    ERBNum(1) = 21.4 * log10(1 + (0.00437 * CUTS(1)));
    ERBNum(2) = 21.4 * log10(1 + (0.00437 * CUTS(2)));
    % Calculate center frequency of channels
    erbCenters = linspace(ERBNum(1),ERBNum(2),CHANNELS);
    % Center frequency by filter number
    fCenter = (10.^(erbCenters/21.4)-1)/4.37e-3;
    % Calculate bandwidths for each ERB band
    for (ff = 1:length(fCenter))
        % Calculate bandwidths
        BW(ff) = 24.7 * (4.37*(fCenter(ff)/1000) + 1);
        fLow(ff) = fCenter(ff) - (BW(ff)/2);
        fUp(ff) = fCenter(ff) + (BW(ff)/2);
    end
    % Remove channels over the Nyquist edge
    Q = fUp < (fs/2);
    if (sum(Q) < CHANNELS)
        warning(['Number of channels reduced: ' num2str(sum(Q))]);
        BW = BW(Q);
        fLow = fLow(Q);
        fCenter = fCenter(Q);
        fUp = fUp(Q);
    end
end



%% FILTERING
%====================================================
% Number of bands
nbands = length(fCenter);

% Nyquist frequency
nyq = fs/2;
len = length(sound);
fsound = zeros(len,length(fCenter));

% BUTTERWORTH
if (strcmp(FILTTYPE,'butter'))
    % FOR each center frequency
    for ff = 1:length(fCenter)
        % Get the Nth-order butterworth filter
        [B,A] = butter(ORDER,[fLow(ff)/nyq,fUp(ff)/nyq]);
        % Filter waveform
        fsound(:,ff) = filter(B,A,sound);
    end
    
    
    % GAMMATONE
elseif (strcmp(FILTTYPE,'gamma'))
    gammaLength = 2^nextpow2(0.128*fs); % At least 128 ms
    timeTmp = (0:gammaLength-1)/fs;
    % Adjust bandwidths, not sure why these are wider for gammatones.
    BW = 1.019.*BW;
    % Prepare arrays
    gammaTone = zeros(gammaLength,length(fCenter));
    timeLead = zeros(length(fCenter));
    phase = 0;
    % Calculate gains, based on integral of impulse
    gains = ((1.019.*BW.*((2*pi)/fs)).^ORDER)./6;
    
    % Calculate impulse response for each channel
    for (ff = 1:length(fCenter))
        if PHASEALIGN
            timeLead(ff) = (ORDER-1)./(2*pi*BW(ff));
            phase = -2*pi*fCenter(ff)*timeLead(ff);
        end
        gammaTone(:,ff) = gains(ff) * fs^3*timeTmp.^(ORDER-1) .* ...
            exp(-2*pi*BW(ff)*timeTmp).*cos(2*pi*fCenter(ff)*timeTmp+phase);
    end
    
    % Filter waveform
    fsound = fftfilt(gammaTone,repmat(sound,1,length(fCenter)));
    
    if (PHASEALIGN)
        % Delay due to time lead
        delay = round(timeLead.*fs);
        % Remove time lead
        for ff = 1:length(fCenter)
            fsound(:,ff) = [fsound(delay(ff)+1:end,ff); zeros(delay(ff),1)];
        end
    end    
end


freqs = fCenter;



%% SUMMED
%====================================================
if (nargout > 2)
    summed = sum(fsound');
    summed = summed';
end



end














