function calcDTF(HRTF);
%function calcDTF(locations, spect, freqs, varargin);
%
% Plots the frequency spectrum across azimuth locations in a colored
%   polar plot.
%
%  INPUTS:
%   locations - Array of azimuth angles (from 0 to 360).
%   spect(f,l) - Spectrum energy with rows represented by 'freqs' bands
%       and each column refers to a different location in the same order 
%       as 'locations'.
%   freqs - Array of center frequencies for each band in spect.
%
%   'hold' - Set to hold plots, useful for subplots.
%
%  REQUIRES: fourierBands, SIGNAL PROCESSING TOOLBOX
%
% @Author:  Eugene Brandewie (Oct 15th, 2018)


[HRTF(:,ss,dd), freqs] = fourierT(HRIR(:,ss,dd),48000, 'N', 1024);

























