function playwav(Y,FS,NBITS);
%function playwav(Y,FS,NBITS);
%
% Plays a wav stimulus via the sound card.
%
% INPUTS:
%   Y = Sampled audio matrix
%   FS = (optional) Sampling frequency [Default: 44,100 Hz]
%   NBITS = (optional) Bit-rate of playback [Default: 16]
%
% @Author:  Eugene Brandewie 6/6/2007
%

if nargin < 3
    NBITS = 16;
end

if nargin < 2
    FS = 44100;
end

p = audioplayer(Y,FS,NBITS);
playblocking(p);
stop(p);


