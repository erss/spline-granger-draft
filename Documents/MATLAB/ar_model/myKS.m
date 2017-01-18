function myKS( y, yhat)
% myKS inputs two signals of equal length and calculates the cumulative
% distribution function of each.  Plots them against each other and
% includes 95% confidence bounds.
%


xbd(1) = min([min(y),min(yhat)]);
xbd(2) = max([max(y),max(yhat)]);

xaxis = xbd(1):0.5:xbd(2);
N1 = length(y);
N2 = length(yhat);
for ii = 1:length(xaxis)
    y_cdf(ii) = length(find(y<xaxis(ii)))/N1;
    yhat_cdf(ii) = length(find(yhat<xaxis(ii)))/N2;
    
end

plot(y_cdf,y_cdf,'r','LineWidth',2);
hold on;
plot(y_cdf,y_cdf-1.36/sqrt(N1),'--r','LineWidth',2);
plot(y_cdf,y_cdf+1.36/sqrt(N1),'--r','LineWidth',2);
plot(y_cdf,yhat_cdf,'k','LineWidth',1.5);
title('KS plot','FontSize',20);
legend('1-1 line', 'Lower CI', 'Upper CI','KS line');
xlim([0 1]);
ylim([0 1]);

% plot(xaxis,y_cdf,'r')
% hold on;
% plot(xaxis,yhat_cdf,'k')


end

