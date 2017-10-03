
model_vals = [30,50,70,90,100,150,400];
figure;
xaxis = -1:.001:3;
nobs = 1000 - model_vals;


for i = 1:length(model_vals)
    
    model_order = model_vals(i);
    nobservations = nobs(i);
   
    subplot 211
    plot(xaxis,fpdf(xaxis,model_order,nobservations-9*model_order),'LineWidth',1.5);
    hold on
    
    subplot 212
    cpts = model_order/3 + 2;
    plot(xaxis,fpdf(xaxis,cpts,nobservations-9*cpts),'LineWidth',1.5);
    nobservations-9*cpts
    hold on
    
end
subplot 211
h=legend('30 lags, fpdf(x,30,700)','50 lags, fpdf(x,50,500)','70 lags, fpdf(x,70,300)','90 lags, fpdf(x,90,100)','100 lags, fpdf(x,100,0)','150 lags, fpdf(x,150,-500)');
set(h,'FontSize',15)
title('Standard Granger')
subplot 212
h=legend('30 lags, fpdf(x,12,862)','50 lags, fpdf(x,18.66,782)','70 lags, fpdf(x,25.33,702)','90 lags, fpdf(x,32,622)','100 lags, fpdf(x,35.33,582)','150 lags, fpdf(x,52,382)');
title('Spline Granger')
set(h,'FontSize',15)


%%
% 
T = 2; 
model_history = 50; 
deltat = 5;
xaxis = -1:.001:3;
fs = 500;

nobs = T*fs - model_history;
cntrlpts = 0:deltat:50;
numcpts = length(cntrlpts) +2;
figure;
plot(xaxis,fpdf(xaxis,numcpts,nobs-64*numcpts),'LineWidth',1.5);
hold on

%%
T = 3; 
model_history = 50; 
deltat = 3;
fs = 500;

nobs = T*fs - model_history;
cntrlpts = [0:deltat:50 50];
numcpts = length(cntrlpts) +2;

plot(xaxis,fpdf(xaxis,numcpts,nobs-64*numcpts),'LineWidth',1.5);

legend('T=2, nobs=950, numspts = 13','T=3, nobs=1450, numspts = 20')



