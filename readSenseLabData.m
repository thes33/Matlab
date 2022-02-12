function [DATA] = readSenseLabData(filename)
% [DATA] = readSenseLabData(filename)
%===============================================================
%
%  Reads the SenseLabOnline data file (csv), ';' deliminated.
%   and returns a DATA structure for further analysis.
%===============================================================
%
% @Author: Eugene Brandewie
% @Date: July 2nd, 2021
%



% Read data file
LINES = readlines(filename);

% Process Header
HEADERS = LINES{1};
LABELS = strsplit(HEADERS,';')';

% Process Data
%==============================================
DATA = struct();

% FOR each data line
for (ll = 2:length(LINES)) % ignore header
    C = strsplit(LINES{ll},';')';
    % IF contains data
    if (~isempty(C{1}))
        
        % FOR each column
        for (hh = 1:length(LABELS))
            
            VAL = C{hh};
            VAL = replace(VAL,',','.'); % Swap from Danish notation
            
            % IF data is string
            if (isempty(str2num(VAL)))
                % Add to cell array
                DATA.(LABELS{hh})(ll-1) = {VAL};
            else
                % Add to numerical array
                DATA.(LABELS{hh})(ll-1) = str2num(VAL);
            end
                
            
        end
        
    end
end









