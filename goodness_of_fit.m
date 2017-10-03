%% Residuals
% DW test
% autocorr
gof_residuals(model_spline);
title('Spline Fit')
gof_residuals(model_standard);
title('Standard Fit')


%% Spectral Test
% plot spectrum
    [m2fit] = grstat(model_true,model_spline);
    [m3fit] = grstat(model_true,model_standard);
    
nelectrodes = size(model_spline.data,1);
for i = 1:nelectrodes
 figure;
 
subplot 211
[faxis, Sxx] =mySpec(model_true.data(i,:),model_true.sampling_frequency,'tapers');
        plot((faxis),10*log(Sxx),'col','k','LineWidth',2);
        xlim([0 model_true.sampling_frequency/4]);
hold on
[faxis, Sxx] =mySpec(model_spline.signal_estimate(i,:),model_true.sampling_frequency,'tapers');
        plot((faxis),10*log(Sxx),'col','r','LineWidth',2);
        xlim([0 model_true.sampling_frequency/4]);
[faxis, Sxx] =mySpec(model_standard.signal_estimate(i,:),model_true.sampling_frequency,'tapers');
        plot((faxis),10*log(Sxx),'col','g','LineWidth',2);
        xlim([0 model_true.sampling_frequency/4]);     
title('Spectrogram','FontSize',20);
ylabel('Power (dB)','FontSize',18)
xlabel('Frequency (Hz)','FontSize',18)
box off


subplot 212


xaxis = m3fit.xaxis(i,:);
plot(xaxis,m3fit.estimate(i,:),'g','LineWidth',1.5)
hold on
plot(xaxis,m3fit.bound1(i,:),'--g',xaxis,m3fit.bound2(i,:),'--g','LineWidth',1)
xaxis = m2fit.xaxis(i,:);
plot(xaxis,m2fit.bound1(i,:),'--r','LineWidth',1)
plot(xaxis,m2fit.bound2(i,:),'--r','LineWidth',1)
plot(xaxis,m2fit.estimate(i,:),'r',xaxis,m2fit.true(i,:),'k','LineWidth',1.5)
axis tight
box off

title('Integrated Spectrum Test','FontSize',20);
ylabel('Cumulative Density','FontSize',18)
xlabel('Averaged Spectrum (1/Hz)','FontSize',18)
% a = get(gca,'XTickLabel');
% set(gca,'XTickLabel',a,'fontsize',16)
axis tight
box off

end



%% Coefficients
% check if the coefficients match
 gof_bootstrap2(model_true,model_spline,model_standard);


