function mySpec( x, f0 )
% MYSPEC computes the spectrum of a signal, x, at a sampling frequency, f0.
%  

dt = 1/f0;
T = length(x)/f0;  % total length of recording (seconds)

df = 1/T;   % frequency resolution
fNQ = f0/2; % Nyquist frequency


xf = fft(x);                  % Fourier transform of x                
Sxx = (2*dt^2/T)*(xf.*conj(xf));  % Power spectrum
Sxx = Sxx(1:length(x)/2+1);   % Remove negative frequencies 

faxis = 0:df:fNQ;             % Frequency axis

plot(faxis,Sxx);     
xlim([0 f0/4]);
xlabel('Frequency (Hz)','FontSize',15);
ylabel('Power','FontSize',15);
title('Spectrogram','FontSize',15);


end

