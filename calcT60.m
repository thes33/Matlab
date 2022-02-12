function [t60] = calcT60(IR, fs, dB_Start, dB_Stop, plot_title);
% function [t60] = calcT60(IR, dB_Start, dB_Stop, fs, plot_title);
%
%  Estimates the Reverberation Time (T60) using the given impulse response
%     and starting and stopping points (dB) for the regression line.


% @Author:
%   Pavel Zahorik 1/25/2008, revised 6/2/2011, revised 2/2/2012
%  Modified - Jun 25, 2020 - Eugene Brandewie
%

PLOT = true;
if (nargin < 5)
    PLOT = false;
end

if nargin < 4
    dB_Start = -5;
    dB_Stop = -35;
end



Decay = flip(10*log10(cumsum(flip(IR.^2))));    % Decay curve (in dB) via backwards integration
Decay = Decay-max(Decay);                           % Set max level to 0 dB

x1 = max(find(Decay > dB_Start));
x2 = max(find(Decay > dB_Stop));
if isempty(x1) | isempty(x2),
    warning('Data exceed dB_Start / dB_Stop range.  T60.m');
end
warning off;
B = regress(Decay(x1:x2),[ones(x2-x1+1,1) [x1:x2]']);
t60 = (((Decay(x1)-60) - B(1))./B(2))./fs;
warning on;

if PLOT
    figure();
    plot([1:length(Decay)]./fs,Decay,'b','LineWidth',[2])
    title(plot_title);
    axis([0 length(IR)/fs -100 0]);
    xlabel('Time (s)');
    ylabel('Energy (dB)');
    set(gca,'ytick',[-100:5:0]);
    grid on
    hold on
    plot([1:length(Decay)]./fs,Decay,'b','LineWidth',[2])
    plot(x1./fs,Decay(x1),'og','LineWidth',[2]);
    plot(x2./fs,Decay(x2),'og','LineWidth',[2]);
    plot([0,5],[B(2)*(0)+B(1), B(2)*(5.*fs)+B(1)],'g','LineWidth',[2])
    text(.3,-10,sprintf('T60 = %5.2f s',t60),'BackgroundColor',[1 1 1]);
    hold off
end