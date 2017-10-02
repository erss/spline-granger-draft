function [faxis, Sxx] = mySpec( x, f0,tapers)
% MYSPEC computes the spectrum of a signal, x, at a sampling frequency, f0.
%  Uses tapers if taper
%  ='tapers', doesn't if taper = 'notapers' (DEFAULT)


% if ~exist('plot_flag','var') % default is to plot
%     taper = 'yesplot';
% end

if ~exist('tapers','var') % default is to not use tapers
    tapers = 'notapers';
end


% if ~exist('color','var') % default is black
%     color = 'k';
% end

if strcmp(tapers,'notapers') %%% -----NO TAPERS------------
    x = x - ones(size(x)).*mean(x);
    dt = 1/f0;
    T = length(x)/f0;  % total length of recording (seconds)
    
    df = 1/T;   % frequency resolution
    fNQ = f0/2; % Nyquist frequency
    
    
    xf = fft(x);                  % Fourier transform of x
    Sxx = (2*dt^2/T)*(xf.*conj(xf));  % Power spectrum
    Sxx = Sxx(1:length(x)/2+1);   % Remove negative frequencies
    
    faxis = 0:df:fNQ;             % Frequency axis
    
    
%     if strcmp(plot_flag,'yesplot') % plots spectrum
%         plot((faxis),10*log(Sxx),'col',color,'LineWidth',2);
%         %         hold on;
%         %         plot(faxis,0.0625./(faxis.^.33))
%         xlim([0 f0/4]);
%         xlabel('Frequency (Hz)','FontSize',16);
%         ylabel('Power (dB)','FontSize',16);
%         title('Spectrogram','FontSize',16);
%     end
%     
    
elseif strcmp(tapers,'tapers') %%% ----USE TAPERS -----------
    %%% TAPER METHOD
    x = x - ones(size(x)).*mean(x);
    dt = 1/f0;
    TW = 4; % 50 ?? time-bandwidth product"
    
    [Sxx, faxis] = pmtm(x,TW,length(x),f0);
    
%     if strcmp(plot_flag,'yesplot') % plots spectrum
%        % plot(faxis,10*log(Sxx),'k','LineWidth',1.5);
%        plot(faxis,10*log(Sxx),'col',color,'LineWidth',1);
% 
%         xlim([0 f0/4]);
%         xlabel('Frequency (Hz)','FontSize',14);
%         ylabel('Power (dB)','FontSize',14);
%         title('Spectrogram','FontSize',14);
%     end
end
end

