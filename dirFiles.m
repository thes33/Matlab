function [FILES] = dirFiles(directory, pattern)
% function dirFiles(directory, pattern)
%
%  Returns a cell array of all the files in the given diretcory that
%  contain the given string pattern.  If no string pattern is given,
%  all files will be returned.
%
%  @Author:  Eugene Brandewie
%  Jun 25th, 2021


% Prepare Inputs
%-------------------------------
USE_PATTERN = true;
if (nargin < 2)
    pattern = [];
    USE_PATTERN = false;
end

% Gather files
%-------------------------------
temp = dir(directory);
cnt = 0;
FILES = {};
for (ff = 3:length(temp))
    fname = {temp(ff).name};
    
    % Pattern
    if (USE_PATTERN)
        if (contains(fname,pattern))
            cnt = cnt +1;
            FILES(cnt) = fname;
        end        
        
    else % No Pattern
        cnt = cnt +1;
        FILES(cnt) = fname;
    end
    
end
% Easier to read output
FILES = FILES';
