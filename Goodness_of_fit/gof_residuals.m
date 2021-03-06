function [notwhite, dw] = gof_residuals( model)
% GOF_RESIDUALS analyzes the goodness of fit of the model residuals. For
% each signal, it plots the true signal and the estimated signal, the model
% residuals over time and an autocorrelation plot of the residuals.
%
% INPUTS:
% . model = structure containing data and estimated data.
%
% OUTPUTS:
% . notwhite = signals that fail the Durbin-Watson test for serial
%              autocorrelation, implying poor model fit.
% . dw       = the Durbin-Watson statistic
data  = model.data;
nelectrodes = size(data,1);
model_order = model.estimated_model_order;
datap = data(:,model_order+1:end);
yestimate = model.signal_estimate;
dt = 1/model.sampling_frequency;
T = model.T;
taxis = dt:dt:T';


residuals = zeros(nelectrodes,length(datap));

for electrode = 1:nelectrodes
   h= figure;
    true_signal =  datap(electrode,:);
    est_signal = yestimate(electrode,:);
    subplot 311
    plot(taxis(model_order+1:end),true_signal,'k');
    %plot(true_signal,'k');
    hold on;
    plot(taxis(model_order+1:end),est_signal,'r');
    %plot(est_signal,'r');
    xlabel('Time(s)');
    legend('true signal','estimated signal')

    residuals(electrode,:) = true_signal-est_signal;
    subplot 312
    plot(taxis(model_order+1:end),residuals(electrode,:),'.');
   plot(residuals(electrode,:),'.');
    title('Residuals')
    xlabel('Time(s)');
    subplot 313
    autocorr(residuals(electrode,:));
    
    suptitle(num2str(electrode))
end

 [dw, pval] = whiteness(data,residuals);

 sig = significance(pval,0.05,'FDR');
notwhite = find(sig);
if isempty(notwhite)
    fprintf('all residuals are white by Durbin-Watson test at significance %g\n',0.05);
else
    fprintf(2,'WARNING: autocorrelated residuals at significance %g for variable(s): %s\n',0.05,num2str(notwhite));
end


end

