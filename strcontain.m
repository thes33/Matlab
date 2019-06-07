function [B, idx] = strcontain(strng, pattern)
% function B = strcontain(strng, pattern)
%
% Returns 'true' (1) if the given pattern exists in the given string.
%  (optional) idx - returns the index of first occurance of the pattern
%     in string (via strfind).  Empty if not found.
%
% @Author: Eugene Brandewie
% @Version: 1.0 - May 28 2018

idx = strfind(strng,pattern);
B =  ~(isempty(idx));