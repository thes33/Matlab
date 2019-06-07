function [b] = hasString(string, pattern)
% function hasString(string, pattern)
%
% Returns 'true' if the given pattern string is present in string.
%


if (~isempty(strfind(string,pattern)))
    b = true;
else b = false;
end