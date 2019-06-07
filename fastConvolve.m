function [output] = fastConvolve(input, impulse);
% function [output] = fastConvolve(input, impulse);
%
%  Fast Fourier Convolution in the frequency domain for use
%      with virtual room impulse responses.
%
% INPUTS:
%   input   - Input data waveform.
%   impulse - Impulse response to convolve data with.
%
% OUTPUT:
%   output  - Output waveform.
%       (Matrix length equal to length of input + length of impulse - 1)
%
%  @Author: Eugene Brandewie
%--------------------------------------------------
%  @Version: 1.0 - (9/19/2013)
%  @Version: 1.1 - (3/28/2018) 
%     - Fixed multi-channel impulse responses.
%     - Added output matrix preallocation.


if (size(input,2) > 1)
    error('Input signal must be columnar.');
end

output = zeros(length(input) + size(impulse,1) -1, size(impulse,2));
for (cc = 1:size(impulse,2))
    
    % Calculate final length
    fLength = length(input) + length(impulse(:,cc)) -1;
    
    % Find next power of 2 from length for FFT
    length2 = pow2(nextpow2(fLength));
    
    % Get Fast-Fourier-Transform of input matrix
    inputF = fft(input, length2);
    
    % Get Fast-Fourier-Transform of impulse matrix
    impulseF = fft(impulse(:,cc), length2);
    
    % Multiply in frequency domain (equals convolution in time domain)
    outputF = inputF .* impulseF;
    
    % Convert back to time domain using inverse FFT
    outputT = real(ifft(outputF, length2));
    
    % Return only the proper matrix length
    output(:,cc) = outputT(1:fLength);
end





