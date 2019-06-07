function wave=nband_gen(low,high,dur,fs)
% function wave=nband_gen(low,high,dur,fs)
%	 wave=nband_gen(low,high,dur,fs)
%   Generate Narrow band Noise with Rayleigh amplitudes and uniform phases

fundemental=1/dur;
low_int=round(low/fundemental);
high_int=round(high/fundemental);
f=low_int:high_int;
len=length(f);
totlen=(dur*fs);
x=randn(1,len);
y=randn(1,len);
spect=[zeros(1,low_int) complex(x,y) zeros(1,totlen-(((2.*high_int)+1))) conj(complex(x(len:-1:1),y(len:-1:1))) zeros(1,low_int-1) ];
length(spect);
wave=real(ifft(spect));

%conj(complex(x,y)) zeros(1,low_int-2)
