function [waveformV, fLow, fUp] = vocoder(waveform, fs, varargin)
% function [waveformV, fLow, fUp] = vocoder(waveform, fs, varargin)
%
%  Channel vocoder processor.
%
%
%   INPUT:
%     waveform = waveform to process, columnar
%     fs   = sampling frequency of waveform
%
%   FILTER-BAND INPUTS:  [Default: 1/6th octave bands]
%     Choose channels based on octaves or fractional octaves:
%       'bandsize', [N/S]  = Bandwidth for each channel filter in octaves. One always centered on 1000 Hz.
%           Variable number of bands based on 'bandwidth' limits.

%     Or choose based on specified band-centers
%       'bands', M = Matrix of band centers (column 1) and bandwidths (column 2). Overrides 'bandsize' option.

%     Or choose equally spaced channels by specifying a number of channels only. ERB-spaced channels.
%       'channels', N = Number of channels spread evenly in 2.1 log-space (within 'bandwidth' limits).
%
%   OPTIONAL INPUTS:
%     'bandwidth', [N1 N2]  = Total bandwidth cut offs [lower higher] in Hz. [Default: 100 8000]
%     'filter', 'S', C = Filter type for band-filters. [Default: 'butter',{3}]: Options: 'butter', 'elliptic'
%          Followed by cell array of filter parameters. Ex: 'butter', {3}; 'elliptic', {6,0.1,80};
%     'envelope', N = Upper frequency cut-off for envelope extraction. [Default: 160 Hz].
%     'carrier', S = Carrier for the vocoded signal: Noise: 'white', 'pink'.  [Default: 'white'].
%
%   OUTPUT:
%       waveformV = Vocode processed output waveform
%
%   OPTIONAL OUTPUTS:
%       fLow = Frequency lowerbounds for channels.
%       fUp = Frequency upperbounds for channels.
%
%  REQUIRES: Signal Processing Toolbox
%
%
% @Author: Eugene Brandewie
% @Version: 1.0 - Jun 5, 2019
%



%% DEFAULT PARAMETERS
%===========================================================
CHANNELS = 8; % Number of channels to vocode.
BANDSIZE = 6; % Size of each individual band, fraction of an octave
BANDWIDTH = [100 8000]; % Bandwidth limits [lower higher] in Hz.
BANDS = []; % Matrix of band centers (column 1) and bandwidths (column 2). Overrides bandsize option.
FILTER_TYPE = 'butter';  % Filter type for band-filters.
FILTER_PARAMS = {3}; % Cell array of filter options.
ENVELOPE = 160; % Upper frequency cut-off for envelope extraction.
CARRIER = 'white'; % Carrier for the vocded signal.

%% INPUT PARAMETERS
%===========================================================

if nargin >= 3
    for ii = 1:nargin-2
        S = varargin{ii};
        if (~isnumeric(S))
            % BANDSIZE
            if strcmpi(S,'bandsize')
                if (~isnumeric(varargin{ii+1}))
                    % OCTAVE
                    if strcmpi(varargin{ii+1},'octave')
                        BANDSIZE = 1;
                    end
                    % 1/2 OCTAVE
                    if strcmpi(varargin{ii+1},'half')
                        BANDSIZE = 2;
                    end
                    % 1/3 OCTAVE
                    if strcmpi(varargin{ii+1},'third')
                        BANDSIZE = 3;
                    end
                    % 1/6 OCTAVE
                    if strcmpi(varargin{ii+1},'sixth')
                        BANDSIZE = 6;
                    end
                else
                    BANDSIZE = varargin{ii+1};
                end
            end
            % BANDS
            if strcmpi(S,'bands')
                if (isnumeric(varargin{ii+1}))
                    BANDS = varargin{ii+1};
                    BANDSIZE = [];
                else
                    error('Bands must be two columns of band centers (column 1) and bandwidths (column 2).');
                end
            end
            % CHANNELS
            if strcmpi(S,'channels')
                if (isnumeric(varargin{ii+1}))
                    CHANNELS = varargin{ii+1};
                    BANDSIZE = [];
                    BANDS = [];
                else
                    error('Bands must be two columns of band centers (column 1) and bandwidths (column 2).');
                end
            end
            
            
            % BANDWIDTH
            if strcmpi(S,'bandwidth')
                BANDWIDTH = varargin{ii+1};
                if (BANDWIDTH(1) < 1 || BANDWIDTH(2) > fs)
                    error('''bandwidth'' must be between 1 and sampling rate.');
                end
            end
            % FILTER TYPE
            if strcmpi(S,'filter')
                FILTER_TYPE = varargin{ii+1};
                if (~iscell(varargin{ii+2}))
                    error('''filter'' must provide cell array of filter input parameters.');
                end
                if (strcmpi(FILTER_TYPE,'butter') || strcmpi(FILTER_TYPE,'elliptic'))
                    FILTER_PARAMS = varargin{ii+2};
                else
                    error('''filter'' must be of type: ''butter'' or ''elliptic''.');
                end
            end
        end
    end
end





%% CREATE BAND-PASS FILTERBANK
%==============================================================

% Nyquist frequency
nyq = fs/2;
% fCenter: list of center frequencies for filters
% fLow: list of low cut-offs for filters
% fUp: list of high cut-offs for filters




% Specified Number of Channels for Bands  [Equally-spaced]
%---------------------------------------------------
if (~isempty(CHANNELS))
    % Channels based on Equivalent Rectangular Bandwidth (ERB) [Moore and Glasberg (1990)]
    
    % Calculate max ERB filter number
    ERBNum(1) = 21.4 * log10(1 + (0.00437 * BANDWIDTH(1)));
    ERBNum(2) = 21.4 * log10(1 + (0.00437 * BANDWIDTH(2)));
    % Calculate center frequency of channels
    erbCenters = linspace(ERBNum(1),ERBNum(2),CHANNELS+1);
    % Center frequency by filter number
    erbFreqs = (10.^(erbCenters/21.4)-1)/4.37e-3;
    % Calculate bandwidths for each ERB band
    disp('Created Channels:');
    for (ff = 2:length(erbFreqs))
        % Calculate bandwidths
        fLow(ff-1) = erbFreqs(ff-1);
        fUp(ff-1) = erbFreqs(ff);
        disp(['   ' num2str(ff-1) ': ' num2str(round(fLow(ff-1))) ' - ' num2str(round(fUp(ff-1))) ' Hz']);
    end
    
    
    % Fractional Octave Bands
    %---------------------------------------------------
elseif (isempty(BANDS))
    if (BANDSIZE == 1) % Octave Band
        fCenter = (10.0^3) * ((2.0) .^ [-10:10]);
        fCenter = fCenter(fCenter>BANDWIDTH(1));
        fCenter = fCenter(fCenter<BANDWIDTH(2));
        fd = (2^(1/2));
        fUp = fCenter .* fd;
        fLow = fCenter ./ fd;
    elseif (BANDSIZE == 2) % 1/2 Octave Band
        fCenter  = (10.0^3) * ((2.0) .^ ([-20:20]./2));
        fCenter = fCenter(fCenter>BANDWIDTH(1));
        fCenter = fCenter(fCenter<BANDWIDTH(2));
        fd = (2^(1/4));
        fUp = fCenter .* fd;
        fLow = fCenter ./ fd;
    elseif (BANDSIZE == 3) % 1/3 Octave Band
        fCenter  = (10.0^3) * ((2.0) .^ ([-30:30]./3));
        fCenter = fCenter(fCenter>BANDWIDTH(1));
        fCenter = fCenter(fCenter<BANDWIDTH(2));
        fd = (2^(1/6));
        fUp = fCenter .* fd;
        fLow = fCenter ./ fd;
    elseif (BANDSIZE == 6)% 1/6th Octave Band
        fCenter  = (10.0^3) * ((2.0) .^ ([-50:40]./6));
        fCenter = fCenter(fCenter>BANDWIDTH(1));
        fCenter = fCenter(fCenter<BANDWIDTH(2));
        fd = (2^(1/12));
        fUp = fCenter .* fd;
        fLow = fCenter ./ fd;
    end
    CHANNELS = length(fCenter);
    
    
    % Custom Bands
    %---------------------------------------------------
else
    CHANNELS = size(BANDS,2);
    % FOR each channel
    for (bb = 1:CHANNELS)
        fCenter(bb)  = BANDS(bb,1);
        fUp(bb) = BANDS(bb,1) - BANDS(bb,2);
        fLow(bb) = BANDS(bb,1) + BANDS(bb,2);
    end
    
    
end %END custom band option





%% CREATE CARRIER SIGNAL
%=====================================================

% White Noise Carrier
%----------------------------------------
if (strcmpi(CARRIER,'white'))
    carrierSignal = randn(length(waveform),1);
    % Normalize RMS
    curRMS = sqrt(mean(carrierSignal.^2));
    carrierSignal = carrierSignal ./ (curRMS/0.5);
    
    % Pink Noise Carrier
    %----------------------------------------
elseif (strcmpi(CARRIER,'pink'))
    % Author: Hristo Zhivomirov,  07/30/13
    % generate white noise sequence
    x = randn(length(waveform),1);
    L = length(x);
    % FFT
    X = fft(x');
    % prepare a vector with frequency indexes
    NumUniquePts = L/2 + 1;     % number of the unique fft points
    k = 1:NumUniquePts;         % vector with frequency indexes
    % manipulate the left half of the spectrum so the PSD
    % is proportional to the frequency by a factor of 1/f,
    % i.e. the amplitudes are proportional to 1/sqrt(f)
    X = X(1:NumUniquePts);
    X = X./sqrt(k);
    % prepare the right half of the spectrum - a conjugate copy of the left
    % one except the DC component and the Nyquist component - they are unique,
    % and reconstruct the whole spectrum
    X = [X conj(X(end-1:-1:2))];
    % IFFT
    p = real(ifft(X));
    % ensure that the length of y is N
    p = p(1, 1:L);
    % form the noise matrix and ensure unity standard
    % deviation and zero mean value (columnwise)
    p = reshape(p, [L, 1]);
    p = bsxfun(@minus, p, mean(p));
    carrierSignal = bsxfun(@rdivide, p, std(p));
    
    
    %TODO: Pulse-train and other carrier stimuli types
end





%% APPLY BAND-PASS FILTERING
%=====================================================

% Preallocate
filtS = zeros(length(waveform),length(fLow));
filtC = zeros(length(carrierSignal),length(fLow));
% FOR each center frequency
for ff = 1:length(fLow)
    
    curOrder = FILTER_PARAMS{1};
    isOrderError = true;
    
    % WHILE filter is failing
    while isOrderError
        
        % Create filter
        if (strcmpi(FILTER_TYPE,'butter'))
            [B,A] = butter(curOrder,[fLow(ff)/nyq,fUp(ff)/nyq]);
            
        elseif (strcmpi(FILTER_TYPE,'elliptic'))
            RP = FILTER_PARAMS{2};
            RS = FILTER_PARAMS{3};
            [B,A] = ellip(curOrder,RP,RS,[fLow(ff)/nyq,fUp(ff)/nyq]);
        end
        isOrderError = false;
        
        % Filter waveform
        filtS(:,ff) = filtfilt(B,A,waveform);
        % Filter carrier
        filtC(:,ff) = filtfilt(B,A,carrierSignal);
        
        % IF failure, reduce filter order
        if ((sum(any(isnan(filtS)))))
            isOrderError = true;
            curOrder = curOrder - 1;
            warning(['Error in filtering band: ' num2str(round(fLow(ff))) ' - ' ...
                num2str(round(fUp(ff))) ': Order reduced to ' num2str(curOrder)]);
        end
        
        
    end
    
    % Calculate waveform RMS for each band
    RMS(ff) = sqrt(mean(filtS(:,ff).^2));
end



%% ENVELOPE EXTRACTION
%=====================================================

% Preallocate
env = zeros(length(waveform),length(fLow));
% FOR each frequency channel
for ff = 1:length(fLow)
    
    % Half-wave rectification
    rec = filtS(:,ff) < 0;
    filtS(rec,ff) = 0;
    
    % Low-pass filtered
    [B,A] = ellip(1,0.1,80,[ENVELOPE/nyq]);
    env(:,ff) = filtfilt(B,A,filtS(:,ff));
    
end




%% MODULATE CARRIER BANDS WITH ENVELOPES
%=====================================================

% Preallocate
modCarrier = zeros(length(waveform),length(fLow));
% FOR each frequency channel
for ff = 1:length(fLow)
    
    % Modulate Band
    modCarrier(:,ff) = filtC(:,ff) .* env(:,ff);
    
    % Apply RMS equalization
    rms = sqrt(mean(modCarrier(:,ff).^2));
    Q = rms / RMS(ff);
    modCarrier(:,ff) = modCarrier(:,ff) ./ Q;
end





%% SUMMATE CARRIER BANDS
%=====================================================

% Preallocate
waveformV = zeros(length(waveform),1);
% FOR each frequency channel
for ff = 1:length(fLow)
    waveformV = waveformV + modCarrier(:,ff);
end










