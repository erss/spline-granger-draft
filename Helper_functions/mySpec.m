function [faxis, Sxx] = mySpec( x, f0, plot_flag )
% MYSPEC computes the spectrum of a signal, x, at a sampling frequency, f0.
%  Plots if PLOT_FLAG = 1, the default

if nargin == 2
    plot_flag = 1;  % Default value 
end
% 
% NO TAPERS
x = x - ones(size(x)).*mean(x);
dt = 1/f0;
T = length(x)/f0;  % total length of recording (seconds)

df = 1/T;   % frequency resolution
fNQ = f0/2; % Nyquist frequency


xf = fft(x);                  % Fourier transform of x                
Sxx = (2*dt^2/T)*(xf.*conj(xf));  % Power spectrum
Sxx = Sxx(1:length(x)/2+1);   % Remove negative frequencies 

faxis = 0:df:fNQ;             % Frequency axis


    if plot_flag == 1
        plot(faxis,Sxx);     
        xlim([0 f0/4]);
        xlabel('Frequency (Hz)','FontSize',15);
        ylabel('Power','FontSize',15);
        title('Spectrogram','FontSize',15);
    end


%%% TAPER METHOD
% x = x - mean(x);
% dt = 1/f0;
% TW = 5; % 50 ??
% 
% [Sxx, faxis] = pmtm(x,TW,length(x),f0);
% 
%     if plot_flag == 1
%         plot(faxis,Sxx);     
%         xlim([0 f0/4]);
%         xlabel('Frequency (Hz)','FontSize',15);
%         ylabel('Power','FontSize',15);
%         title('Spectrogram','FontSize',15);
%     end

end