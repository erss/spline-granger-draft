function [faxis, S ] = myTheoreticalSpectrum( model_coefficients, noise_process, f0)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

T = (length(noise_process))/f0;  % total length of recording (seconds)
df = 1/T;                          % frequency resolution
fNQ = f0/2;                        % Nyquist frequency
faxis = 0:df:fNQ;                  % Frequency axis
dt = 1/f0;
model_order = length(model_coefficients);

noise_variance = var(noise_process);
    for f = 2:length(faxis)
       freq = faxis(f);
       temp  = sum(model_coefficients.*exp(-2*pi*1i*freq.*(1:model_order).*dt));
       temp = abs(1-temp).^2;
       S(f) = noise_variance*dt/temp;
    end
end

