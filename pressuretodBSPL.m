function [dBSPL] = pressuretodBSPL(pressure);
%function [intensity] = dbSPLtoPressure(dbSPL);
%
%  Calculates the sound pressure level in decibels (dB SPL) based on the given pressure (in Pascals).
%
%  @Author: Eugene Brandewie 

dBSPL = 20.*log10(pressure./0.000020);