function [NB] = createNoiseBurst(fs, dur, varargin);
% [NB] = createNoiseBurst(fs, dur, varargin);
% Generates a noise burst using the following parameters.
%
%   INPUTS:
%       fs - Sampling frequency
%       dur - Stimulus duration (seconds).
%
%   OPTIONAL INPUTS:
%       'max', N (num) - Maximum peak value.
%       'rms', N (num) - Output root-mean-square level. (Max peaks will be truncated).
%       'type', S (string) - Spectral energy distribution. ['w'/'white' or 'p'/'pink']
%       'band', [Order Lower Upper] - Uses an [Order] Butterworth
%           filter using the given upper and lower cutoffs. Set Cutoff to
%           -1 to ignore cutoff, ex: [3 -1 3000] for a 3rd order low-pass
%           filtered noise at 3 kHz.
%       'ramp', N (num) - Ramp the noise on/off for the given msec duration with cos2 ramps.
%           Ramps applied after RMS, so final RMS is likely under due to amplitude ramping.
%
%   OUTPUT:
%       NB - Noise burst waveform.
%
% @Author: Eugene Brandewie, Hristo Zhivomirov

% Changelog:
%   Dec 19, 2018 - Added pink noise support.
%   Mar 15, 2021 - Added cos2 ramping option.
%   Nov 10, 2021 - Better messaging concerning RMS and max.


% INPUT PARAMETERS
%==========================================


%% Handle optional input parameters
%-----------------------------------------------------------
RMS = false;
MAX = false;
TYPE = 'w';
BANDPASS = false;
RAMP = 0;

if nargin >= 3
    for ii = 1:nargin-2
        S = varargin{ii};
        if (~isnumeric(S))
            % MAX VALUE
            if strcmpi(S,'max')
                maxValue = varargin{ii+1};
                MAX = true;
            end
            % RMS VALUE
            if strcmpi(S,'rms')
                rmsValue = varargin{ii+1};
                RMS = true;
            end
            % NOISE TYPE
            if strcmpi(S,'type')
                TYPE = varargin{ii+1};
            end
            % BANDPASS FILTER
            if strcmpi(S,'band')
                band = varargin{ii+1};
                order = band(1);
                lower = band(2);
                upper = band(3);
                BANDPASS = true;
            end
            % RAMP DURATION
            if strcmpi(S,'ramp')
                RAMP = varargin{ii+1};
            end
        end
    end
end




% NOISE GENERATION
%==========================================


% White Noise
%--------------------------------
if (strcmp(TYPE,'w') || strcmp(TYPE,'white'))
    NB = randn(round(dur*fs),1);
    
    
    % Pink Noise
    %--------------------------------
elseif (strcmp(TYPE,'p') || strcmp(TYPE,'pink'))
    % Author: Hristo Zhivomirov,  07/30/13
    
    % generate white noise sequence
    x = randn(round(dur*fs),1);
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
    p = bsxfun(@rdivide, p, std(p));
    NB = p;
    
else % ELSE type not found
    error(['Noise generation of type ''' TYPE ''' not found.']);
end


% BANDPASS FILTERING
%============================================
if (BANDPASS)
    % IF Low-pass filter
    if (lower == -1 && upper > 0)
        [B,A] = butter(order,upper/(fs/2),'low');
        NB = filter(B,A,NB);
        
        % IF High-pass filter
    elseif (lower > 0 && upper == -1)
        [B,A] = butter(order,lower/(fs/2),'high');
        NB = filter(B,A,NB);
        
        % IF Band-pass filter
    elseif (lower > 0 && upper > 0)
        [B,A] = butter(order,[lower/(fs/2) upper/(fs/2)],'bandpass');
        NB = filter(B,A,NB);
    end
end


% APPLY MAX VALUE LIMITS
%============================================
if (MAX && ~RMS)
    if (max(abs(NB)) > maxValue)
        preRMS = sqrt(mean(NB.^2));;
        Q = max(abs(NB)) / maxValue;
        NB = NB / Q;     
        postRMS = sqrt(mean(NB.^2));;  
        disp(['   RMS level [' num2str(preRMS) '] adjusted [' num2str(postRMS) '] to limit peaks.']); 
    end
end


% APPLY RMS VALUE LIMITS
%============================================
if (RMS)
    rmsNB = sqrt(mean(NB.^2));
    Q = rmsNB / rmsValue;
    NB = NB / Q;
    
    if (MAX) % Truncate maximums
        % Pos
        idx = find(NB >= maxValue);
        posNum = length(idx);
        NB(idx) = maxValue;
        % Neg
        idx = find(NB <= -maxValue);
        negNum = length(idx);
        NB(idx) = -maxValue;
        
        truncSamps = posNum + negNum;
        if (truncSamps > 0) disp(['   ' num2str(truncSamps) ' / ' num2str(length(NB)) ' samples truncated to max value.']); end
    end
end


% APPLY COS2 RAMPS
%============================================
if (RAMP > 0)
    % Get sample lengths
    sampRamp = round(RAMP*(fs/1000));
    % Create the ramps
    onsetRamp = ((sin(pi*(3/2):pi/(sampRamp-1):pi*(5/2))+1)/2).^2;
    offsetRamp = fliplr((sin(pi*(3/2):pi/(sampRamp-1):pi*(5/2))+1)/2).^2;
    % Create the full-scale middle portion
    middle = ones(1, length(NB)-(length(onsetRamp)+length(offsetRamp)));
    % Create final gain curve
    gains = [onsetRamp, middle, offsetRamp];
    % Apply the gain
    NB = gains' .* NB;
    
end














