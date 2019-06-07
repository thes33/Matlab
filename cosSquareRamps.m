function [outWave] = cosSquareRamps(waveform, fs, onset, offset)
%[outWave] = cosSquareRamps(waveform, fs, onset, offset)
%
%  Applies onset and offset raised cosine gates to a sound. If only
%   a single duration value indicated, applies to both onset and offset.
%
%  INPUTS:
%   waveform - Waveform to modify (must be columnar).
%   fs - Sampling frequency.
%   onset - Onset ramp duration in ms [0 = no ramp].
%   offset (optional) - Offset ramp duration in ms [0 = no ramp].
%
%  OUTPUTS:
%   outWave = modified waveform
%
% @Author: Eugene Brandewie
%
%

%% INPUT PARAMETERS
%==============================================

if (size(waveform,2) > size(waveform,1))
    error('Waveform must be columnar: apply transposition.');
end

if (nargin < 3)
     error('Not enoguh input arguements:  need onset/offset duration.');
elseif (nargin == 3)
    offset = onset;
end

if (onset == 0 && offset == 0)
     error('Either the onset or offset ramp must be greater than 0 ms.');
end


%% INPUT PARAMETERS
%==============================================

% Deal with single samples
if (rem(length(waveform), 2) ~=0 && size(waveform, 2) == 1)
    waveform = waveform(1:end-1);
elseif (rem(length(waveform), 2) ~=0 && size(waveform, 2) == 2)
    waveform = waveform(1:end-1, 1:2);
end

% Get sample lengths
sOnset = round(onset*(fs/1000));
sOffset = round(offset*(fs/1000));
if (sOnset+sOffset) > length(waveform)
    error ('The onset/offset ramps cannot be greater than the waveform.');
end


%% APPLY GAIN RAMPING
%==============================================

% Create the ramps
onsetRamp = ((sin(pi*(3/2):pi/(sOnset-1):pi*(5/2))+1)/2).^2;
offsetRamp = fliplr((sin(pi*(3/2):pi/(sOffset-1):pi*(5/2))+1)/2).^2;
% Create the full-scale middle portion
middle = ones(1, length(waveform)-(length(onsetRamp)+length(offsetRamp)));
% Create final gain curve
gains = [onsetRamp, middle, offsetRamp];

% Apply the gain
%------------------------
% IF 2 channel
if size(waveform, 2)==2
    outWave = [gains' .* waveform(:, 1), gains' .* waveform(:, 2)];
% ELSE 1 channel
else
    outWave = gains' .* waveform;
end









