function [sData] = smoothData(data, varargin)
% function [sData] = smoothData(data, parameters...);
%
% Smooths the given data using a specified filter.
%
% INPUTS:
%   data = Data function to be smoothed, can have multiple columns.
%
% FILTER TYPE:  [Default: 'avg',5]
%   'avg',N = Moving average filter with the given length.
%   'hann',N = Hann window filter with the length.
%
%


%% INPUT PARAMETERS
%=======================================================
if (nargin < 1)
    error('Not enough input arguments.');
end
if (size(data,2) > size(data,1))
    error('Data must be columnar.');
end

HANN = 0;
AVG = 0;
N = 5;

% FOR each parameter
if nargin > 1
    for ii = 1:(nargin-1)
        Q = varargin{ii};
        if (~isnumeric(Q))
            % Hann
            if (strcmp(Q,'hann'))
                HANN = 1;
                N = varargin{ii+1};                
            end
            % Avg
            if (strcmp(Q,'avg'))
                AVG = 1;
                N = varargin{ii+1};                
            end
        end
    end
end

% Defaults
if (nargin < 2)
    AVG = 1; 
    N = 5;
end




%% FILTER COEFICIENTS
%=======================================================

if (AVG)
    B = ones(1,N) * (1/N);
    A = 1;
elseif (HANN)
    B = hann(N) * (1/N);
    A = 1;
end



%% FILTER FUNCTION
%=======================================================

sData = filter(B,A,data);



    
    
    
    
    
    
    
    
    

end % END function

