function [t60, t60m, fbands] = calcT60Bands(IR, fs, varargin)
% function [[t60, t60m, fbands] = calcT60Bands(IR, fs, varargin)
%
%  Estimates the Reverberation Time (T60) using the given impulse response
%     for each frequency band and starting and stopping points (dB) for the regression line.
%
% INPUTS:
%      IR = Impulse response waveform (samples x 1)
%      fs = sampling frequency
%
% OPTIONAL INPUTS:
%      'regress' =  Use the regression line method on decay curve. [Default]
%      'snoise' - Stopped noise method, averaged over 10 samples.
%      'plot' = Plot results by frequency band.
%      'dBStart' = Point on decay curve to start estimate analysis in dB. [Default -5 dB]
%      'dBStop' = Point on decay curve to stop estimate analysis, in dB. [Default -35 dB]
%
% OUTPUTS:
%      t60 = Reverberation Time (T60) for each octave band [125 250 500 1000 2000 4000 8000 16000].
%      t60m = Broadband reverbation time (T60) as average of 500 and 1000 Hz bands.
%      fbands = Frequency band centers: [125 250 500 1000 2000 4000 8000 16000] in Hz.
%
%   REQUIRES SIGNAL PROCESSING TOOLBOX
%
% @Author:
%   Eugene Brandewie - Jul 1, 2020
%   Regression method by Pavel Zahorik Jan 25, 2008
%   Stopped noise method by James M. Kates, Jun 15, 2020.
%


%% Input Handling
%-----------------------------------------------------------



% Regress requirement
if ~(exist('regress.m','file') == 2)
    error('Requires ''regress.m'' for regression analysis.');
end


if (nargin < 1)
    error('Requires an impulse response waveform to analyze.');
elseif (nargin < 2)
    error('Provide the sampling rate of the impulse response.');
end

if (size(IR,2) > 1)
    error('Too many columns. Can only analyze one waveform channel.');
end



%% Handle optional input parameters
%-----------------------------------------------------------

% Set defaults
PLOT = false;
METHOD = 'regress';
dBStart = -5;
dBStop = -35;

if nargin >= 3
    for ii = 1:nargin-2
        S = varargin{ii};
        if (~isnumeric(S))
            % PLOT
            if strcmpi(S,'plot')
                PLOT = true;
            end
            % STOPPED NOISE METHOD
            if strcmpi(S,'snoise')
                METHOD = 'snoise';
            end
            % REGRESSION METHOD
            if strcmpi(S,'regress')
                METHOD = 'regress';
            end
            % ANALYSIS START
            if strcmpi(S,'dBStart')
                dBStart = varargin{ii+1};
            end
            % ANALYSIS STOP
            if strcmpi(S,'dBStop')
                dBStop = varargin{ii+1};
            end
        end
    end
end


%% Create Octave band filters
% Octave Band
CUTS = [124 16001];
fbands = (10.0^3) * ((2.0) .^ [-10:10]);
fbands = fbands(fbands>CUTS(1));
fbands = fbands(fbands<CUTS(2));
fd = (2^(1/2));
fUp = fbands .* fd;
fLow = fbands ./ fd;

% Nyquist frequency
nyq = fs/2;
len = length(IR);
fsound = zeros(len,length(fbands));
ORDER = 3;  % 3rd order butterworth filters


if PLOT
    figure();
    annotation('textbox', [0 0.9 1 0.1], ...
        'String', 'Reverberation Time by Octave Bands', ...
        'EdgeColor', 'none', ...
        'HorizontalAlignment', 'center', ...
        'FontSize',14);
end

% FOR each center frequency
for ff = 1:length(fbands)
    % Get the Nth-order butterworth filter
    [B,A] = butter(ORDER,[fLow(ff)/nyq,fUp(ff)/nyq]);
    % Filter waveform
    fsound(:,ff) = filter(B,A,IR);
    
    
    % Regression Method
    %-------------------------------------------
    if (strcmp(METHOD,'regress'))
        % Calculate decay function
        Decay = flip(10*log10(cumsum(flip(fsound(:,ff).^2))));    % Decay curve (in dB) via backwards integration
        Decay = Decay-max(Decay);                           % Set max level to 0 dB
        
        fitStart = max(find(Decay > dBStart));
        fitStop = max(find(Decay > dBStop));
        if isempty(fitStart) | isempty(fitStop),
            warning('Data exceed dB_Start / dB_Stop range.  T60.m');
        end
        fit = regress(Decay(fitStart:fitStop),[ones(fitStop-fitStart+1,1) [fitStart:fitStop]']);
        t60(ff) = (((Decay(fitStart)-60) - fit(1))./fit(2))./fs;
        
        % PLOTTING
        if (PLOT)
            subplot(4,2,ff);
            plot([1:length(Decay)]./fs,Decay,'b','LineWidth',[2])
            axis([0 length(IR)/fs -60 0]);
            xlabel('Time (s)');
            ylabel('Energy (dB)');
            set(gca,'ytick',[-100:10:0]);
            grid on
            hold on
            plot([1:length(Decay)]./fs,Decay,'b','LineWidth',[2])
            plot(fitStart./fs,Decay(fitStart),'og','LineWidth',[2]);
            plot(fitStop./fs,Decay(fitStop),'og','LineWidth',[2]);
            plot([0,5],[B(2)*(0)+B(1), B(2)*(5.*fs)+B(1)],'g','LineWidth',[2])
            text(.3,-10,sprintf('T60 = %5.3f s',t60(ff)),'BackgroundColor',[1 1 1]);
            hold off
            title([num2str(fbands(ff)) ' Hz']);
        end
        
        % Stopped Noise Method
        %-------------------------------------------
    elseif (strcmp(METHOD,'snoise'))
        
        filtIR = fsound(:,ff);
        fitStart = zeros(10,1);
        fitStop = zeros(10,1);
        clear envs fits;
        
        % FOR 10 iterations
        for (ii = 1:10)
            
            % Create 2-sec Gaussian white noise signal
            rng('shuffle');
            noiseLength = (2*fs);
            rdNoise = randn(noiseLength,1);
            
            % Measure peak
            peakPnt = find(abs(filtIR) == max(abs(filtIR)));
            % Convolve noise with filtered room response
            convNoise = conv(rdNoise,filtIR);
            % Acquire tail
            noiseTail = convNoise(noiseLength + peakPnt : end);
            % Calculate envelope of decay tail in dB
            T = bb_HilbertFIR(noiseTail,fs);
            tailEnv = sqrt((real(T).^2) + imag(T).^2);
            tailEnv = 20*log10(tailEnv/max(tailEnv));
            samps = (1:length(tailEnv))';
            % Polynomial fit parameters, 1st-order
            polyN = polyfit(samps,tailEnv,1);
            % Fitted linear regression curve
            fit = polyN(2) + (polyN(1) * samps);
            % Estimate time points on linear fit
            fitStarttemp = find(fit <= dBStart); 
            fitStarttemp = fitStarttemp(1);
            fitStart(ii) = fitStarttemp;
            % Ensure enough decay for stop point
            fitStoptemp = find(fit <= dBStop);
            if (isempty(fitStoptemp))
                warning(['Data exceeds dB_Start / dB_Stop range. Using ' num2str(fit(end)) ' for dBStop.']);
                fitStop(ii) = length(fit);
                % Modified range factor
                rngFactor = -(60 / (fit(end) - fit(fitStart(ii))));
            else
                fitStoptemp = fitStoptemp(1);
                fitStop(ii) = fitStoptemp;
                % Calculate range factor
                rngFactor = -(60 / (fit(fitStop(ii)) - fit(fitStart(ii))));
            end
            
            % Estimate time to -60 dB based on polynomial within measured range
            T60est(ii) = rngFactor * ((fitStop(ii) - fitStart(ii)) / fs);
            
            % Grab measures to compute averages for graphs
            if (PLOT)
                envs(:,ii) = tailEnv;
                fits(:,ii) = fit;
            end
            
        end
        
        % Calculate mean across samples
        t60(ff) = mean(T60est);
        mStart = mean(fitStart')';
        mStop = mean(fitStop')';
        
        
        % PLOTTING
        if (PLOT)
            % Mean decay function
            decayEnv = mean(envs')';
            mFit = mean(fits')';
            
            subplot(4,2,ff);
            plot([1:length(decayEnv)]./fs,decayEnv,'b','LineWidth',[1])
            axis([0 length(IR)/fs -60 0]);
            xlabel('Time (s)');
            ylabel('Energy (dB)');
            set(gca,'ytick',[-100:10:0]);
            grid on
            hold on
            plot([1:length(mFit)]./fs,mFit,'g','LineWidth',[2])
            plot(mStart./fs,dBStart,'og','LineWidth',[2]);
            plot(mStop./fs,dBStop,'og','LineWidth',[2]);
            plot(t60(ff),-60,'rx','LineWidth',[2]);
            text(.3,-10,sprintf('T60 = %5.3f s',t60(ff)),'BackgroundColor',[1 1 1]);
            hold off
            title([num2str(fbands(ff)) ' Hz']);
            
        end
        
    end % END for method switch
    
    
end % END for each frequency band

% Broadband T60 (mean 500 and 1000 Hz)
t60m = mean([t60(3) t60(4)]);

if (PLOT)
    annotation('textbox', [0 0.87 1 0.1], ...
        'String', ['Broadband t60 = ' num2str(t60m)], ...
        'EdgeColor', 'none', ...
        'HorizontalAlignment', 'center', ...
        'FontSize',12);
end




end








%% Imbedded modified Octave regress function
%---------------------------------

% Copyright (C) 2005, 2006 William Poetra Yoga Hadisoeseno
% Copyright (C) 2011 Nir Krakauer
%
% This program is free software; you can redistribute it and/or modify it under
% the terms of the GNU General Public License as published by the Free Software
% Foundation; either version 3 of the License, or (at your option) any later
% version.
%
% This program is distributed in the hope that it will be useful, but WITHOUT
% ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
% FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
% details.
%
% You should have received a copy of the GNU General Public License along with
% this program; if not, see <http://www.gnu.org/licenses/>.

% -*- texinfo -*-
% @deftypefn {Function File} {[@var{b}, @var{bint}, @var{r}, @var{rint}, @var{stats}] =} regress (@var{y}, @var{X}, [@var{alpha}])
% Multiple Linear Regression using Least Squares Fit of @var{y} on @var{X}
% with the model @code{y = X * beta + e}.
%
% Here,
%
% @itemize
% @item
% @code{y} is a column vector of observed values
% @item
% @code{X} is a matrix of regressors, with the first column filled with
% the constant value 1
% @item
% @code{beta} is a column vector of regression parameters
% @item
% @code{e} is a column vector of random errors
% @end itemize
%
% Arguments are
%
% @itemize
% @item
% @var{y} is the @code{y} in the model
% @item
% @var{X} is the @code{X} in the model
% @item
% @var{alpha} is the significance level used to calculate the confidence
% intervals @var{bint} and @var{rint} (see `Return values' below). If not
% specified, ALPHA defaults to 0.05
% @end itemize
%
% Return values are
%
% @itemize
% @item
% @var{b} is the @code{beta} in the model
% @item
% @var{bint} is the confidence interval for @var{b}
% @item
% @var{r} is a column vector of residuals
% @item
% @var{rint} is the confidence interval for @var{r}
% @item
% @var{stats} is a row vector containing:
%
%   @itemize
%   @item The R^2 statistic
%   @item The F statistic
%   @item The p value for the full model
%   @item The estimated error variance
%   @end itemize
% @end itemize
%
% @var{r} and @var{rint} can be passed to @code{rcoplot} to visualize
% the residual intervals and identify outliers.
%
% NaN values in @var{y} and @var{X} are removed before calculation begins.
%
% @end deftypefn

% References:
% - Matlab 7.0 documentation (pdf)
% - 《大学数学实验》 姜启源 等 (textbook)
% - http://www.netnam.vn/unescocourse/statistics/12_5.htm
% - wsolve.m in octave-forge
% - http://www.stanford.edu/class/ee263/ls_ln_matlab.pdf

function [b, bint, r, rint, stats] = regress (y, X, alpha)

if (nargin < 2 || nargin > 3)
    print_usage;
end

if (~ ismatrix (y))
    error ("regress: y must be a numeric matrix");
end
if (~ ismatrix (X))
    error ("regress: X must be a numeric matrix");
end

if (size(y,2) ~= 1)
    error ("regress: y must be a column vector");
end

if (size(y,1) ~= size(X,1))
    error ("regress: y and X must contain the same number of rows");
end

if (nargin < 3)
    alpha = 0.05;
elseif (~ isscalar (alpha))
    error ("regress: alpha must be a scalar value")
end

notnans = ~ logical (sum (isnan ([y X]), 2));
y = y(notnans);
X = X(notnans,:);

[Xq Xr] = qr (X, 0);
pinv_X = Xr \ Xq';

b = pinv_X * y;

if (nargout > 1)
    
    n = rows (X);
    p = columns (X);
    dof = n - p;
    t_alpha_2 = tinv (alpha / 2, dof);
    
    r = y - X * b; % added -- Nir
    SSE = sum (r .^ 2);
    v = SSE / dof;
    
    % c = diag(inv (X' * X)) using (economy) QR decomposition
    % which means that we only have to use Xr
    c = diag (inv (Xr' * Xr));
    
    db = t_alpha_2 * sqrt (v * c);
    
    bint = [b + db, b - db];
    
end

if (nargout > 3)
    
    dof1 = n - p - 1;
    h = sum(X.*pinv_X', 2); %added -- Nir (same as diag(X*pinv_X), without doing the matrix multiply)
    
    % From Matlab's documentation on Multiple Linear Regression,
    %   sigmaihat2 = norm (r) ^ 2 / dof1 - r .^ 2 / (dof1 * (1 - h));
    %   dr = -tinv (1 - alpha / 2, dof) * sqrt (sigmaihat2 .* (1 - h));
    % Substitute
    %   norm (r) ^ 2 == sum (r .^ 2) == SSE
    %   -tinv (1 - alpha / 2, dof) == tinv (alpha / 2, dof) == t_alpha_2
    % We get
    %   sigmaihat2 = (SSE - r .^ 2 / (1 - h)) / dof1;
    %   dr = t_alpha_2 * sqrt (sigmaihat2 .* (1 - h));
    % Combine, we get
    %   dr = t_alpha_2 * sqrt ((SSE * (1 - h) - (r .^ 2)) / dof1);
    
    dr = t_alpha_2 * sqrt ((SSE * (1 - h) - (r .^ 2)) / dof1);
    
    rint = [r + dr, r - dr];
    
end

if (nargout > 4)
    
    R2 = 1 - SSE / sum ((y - mean (y)) .^ 2);
    %    F = (R2 / (p - 1)) / ((1 - R2) / dof);
    F = dof / (p - 1) / (1 / R2 - 1);
    pval = 1 - fcdf (F, p - 1, dof);
    
    stats = [R2 F pval v];
    
end

end




%% Imbedded Modified Hilbert function
%---------------------------------

function xenv = bb_HilbertFIR(x,fsamp)
% Function to extract the speech envelope using a time-domain
% implementation of the Hilbert transform. The MATLAB hilbert function
% uses a FFT-based algorithm which can generate artifacts when there is a
% gap between sentences. The implementation here designs a linear-phase
% 90-deg phase shift filter which is used for a time-domain HT.
%
% Caling arguments:
% x       signal matrix [nsamp,nband]
% fsamp   sampling rate in Hz
%
% Returned values:
% xenv    envelope of the signal in each band [nsamp,nband]
%
% James M. Kates, 4 May 2020.

% Input speech parameters
nsamp=size(x,1);
nband=size(x,2);

% Design the Hilbert transform filter
tHil=16; %Filter length in msec
nHil=0.001*tHil*fsamp; %Filter length in samples
nHil=2*round(nHil/2); %Force even filter length
nHil2=nHil/2; %Transient duration for the linear-phase filter
fNyq=0.5*fsamp; %Nyquist rate
fLow=60; %Filter lower edge in Hz
h=firpm(nHil,[fLow/fNyq,0.99],[1,1],'Hilbert'); %Design the FIR filter

% Loop to compute the envelope in each band
xenv=zeros(size(x));
for k=1:nband
    xreal=x(:,k); %Extract the signal in the band
    ximag=conv(xreal,h); %Imaginary part of the analytic signal
    ximag=ximag(nHil2+1:nHil2+nsamp); %Remove the transients
    env=sqrt(xreal.^2 + ximag.^2); %Envelope magnitude from HT
    env=max(env,eps); %Ensure no zero or negative values
    xenv(:,k)=env; %Save the envelope
end
end













