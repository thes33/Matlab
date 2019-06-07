function [pressure] = dBSPLtoPressure(dbSPL);
%function [intensity] = dbSPLtoPressure(dbSPL);
%
%  Calculates the pressure (in Pascals) based on the given sound pressure level in decibels (dB SPL).
%
%  @Author: Eugene Brandewie 

pressure = 0.000020 .* (10.^(dbSPL./20));