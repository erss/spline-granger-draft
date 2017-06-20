%%%% goodness_of_fit_residuals

datap = data(:,model_order+1:end);
%datap = data(:,1:end-model_order);

residuals = zeros(nelectrodes,length(datap));

for electrode = 1:nelectrodes
    figure;
    true_signal =  datap(electrode,:);
    est_signal = yestimate(electrode,:);
    subplot 311
    %plot(taxis(model_order+1:end),true_signal,'k');
    plot(true_signal,'k');
    hold on;
    %plot(taxis(model_order+1:end),est_signal,'r');
    plot(est_signal,'r');
    xlabel('Time(s)');
    legend('true signal','estimated signal')

    residuals(electrode,:) = true_signal-est_signal;
    subplot 312
   % plot(taxis(model_order+1:end),residuals(electrode,:),'.');
   plot(residuals(electrode,:),'.');
    title('Residuals')
    xlabel('Time(s)');
    subplot 313
    autocorr(residuals(electrode,:));
    
    suptitle(num2str(electrode))
end

[dw pval] = whiteness(data,residuals)
% A standard rule of thumb is that |dw < 1| or |dw > 3| indicates a high
% chance of residuals serial correlation; this implies poor VAR model fit.
sig = significance(pval,0.05,'FDR')
