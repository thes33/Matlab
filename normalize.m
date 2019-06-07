function [new] = normalize(waveform, N);
% Normalizes a waveform on a N to -N scale (default: 1.0);

if (nargin < 2)
	N = 1.0;
end

new=waveform/(N*max(max(abs(waveform))));
