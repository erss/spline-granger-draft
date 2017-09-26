%%%%% single node analysis
%%% set model coefficitens to single_node_order20 & single_node_low_freq
clear all;
ntrials = 1000;
ct_spline = zeros(1,ntrials);
ct_standard = zeros(1,ntrials);
ts_spline = zeros(1,ntrials);
ts_stand = zeros(1,ntrials);
dwstandard = zeros(1,ntrials);
dwspline = zeros(1,ntrials);
fails =[];
fails_standard =[];

   
    
    
%%% Model type ------------------------------------------------------------
model_true.noise_type = 'white'; % 'white', 'pink', 'real'

%%% Simulation parameters -------------------------------------------------

model_true.sampling_frequency = 500;
model_true.T = 2;   % time in seconds of window
model_true.noise = 0.25;
taxis = (1/model_true.sampling_frequency):(1/model_true.sampling_frequency):model_true.T;
model_true.taxis = taxis;

    %load('large_network_coef.mat');
   % load('ninenode_exp_stand.mat');
    model_true.true_coefficients = single_node_low_freq; %%%% MODIFY COEFFICIENTS HERE!
    model_true.model_coefficients = model_true.true_coefficients;   

%%% Define model inputs for spline Granger & standard Granger -------------

model_true.s = 0.5;                     % tension parameter for spline
model_true.estimated_model_order = 30;  % model_order used to estimate

number_of_knots      = floor(model_true.estimated_model_order/3);
model_true.cntrl_pts = [0:5:30]; %make_knots(model_true.estimated_model_order,number_of_knots);

%%% Define network testing parameters -------------------------------------

model_true.q = 0.05;            % FDR max number acceptable proportion of false discoveries
model_true.nsurrogates = 1000;   % number of surrogates used for bootstrapping
model_true.nrealizations = 20; % number of realizations used for spectral testing



for i = 1:ntrials

    simulate_network;
    infer_network;
    
    ct_spline(i)  = model_spline.computation_time;
    ct_standard(i) = model_standard.computation_time;
    
    
    acc_spline(i)  = model_spline.accuracy;
    acc_standard(i) = model_standard.accuracy;
    
    [nw, dwstandard(i),~] = dwstat(model_standard);
    fails_standard = [fails_standard nw];
    [nw, dwspline(i),~] = dwstat(model_spline);
    fails = [fails nw];
    
    [m2fit, m3fit] = grstat1(model_true,model_spline,model_standard);
 ts_spline(i) = m2fit.stat;
 ts_stand(i)  = m3fit.stat;
end

%%

%%% Fig 1
%%% Plot signal trace
subplot(3,2,1)
plot(model_true.taxis,model_true.data(1,:),'k', 'LineWidth', 1);
hold on;
plot([0,.2], [min(model_true.data(1,:)), min(model_true.data(1,:))], 'k', 'LineWidth', 2.5);
box off
set(gca,'YTickLabel',[]);
set(gca,'XTickLabel',[]);
set(gca,'XTick',[]);
set(gca,'YTick',[]);
title('Simulated Signal','FontSize',20);
xlabel('Time (s)','FontSize',18);

%%% Plot spectrogram
subplot(3,2,2)
mySpec(model_true.data(1,:),model_true.sampling_frequency,'yesplot','tapers');
a = get(gca,'YTickLabel');
set(gca,'YTickLabel',a,'fontsize',16)
title('Spectrogram','FontSize',20);
ylabel('Power (dB)','FontSize',18)
xlabel('Frequency (Hz)','FontSize',18)
box off
axis tight

%%% Plot coefficients
subplot(3,2,[3 4])
gof_bootstrap(model_true,model_spline,model_standard);
box off
title('AR Coefficients','FontSize',20);
ylabel('Magnitude','FontSize',18)
xlabel('Lag (s)','FontSize',18)
% a = get(gca,'XTickLabel');
% set(gca,'XTickLabel',a,'fontsize',16)
plot([0 0.06],[0 0],'color',[.57 .57 .57],'LineWidth',1.7)
        h = legend('Spline Estimated','Standard Estimated','True ');
        set(h,'FontSize',15,'Location','NorthEast')

%%% Plot spectral test
subplot(3,2,5)
%gof_spectrum(model_true,model_spline,model_standard);
xaxis = m3fit.xaxis;
plot(xaxis,m3fit.estimate,'g','LineWidth',1.5)
hold on
plot(xaxis,m3fit.bound1,'--g',xaxis,m3fit.bound2,'--g','LineWidth',1)
xaxis = m2fit.xaxis;
plot(xaxis,m2fit.bound1,'--r',xaxis,m2fit.bound2,'--r','LineWidth',1)
plot(xaxis,m2fit.estimate,'r',xaxis,m2fit.true,'k','LineWidth',1.5)
% 
% 
% set(gca,'YTickLabel',[]);
% set(gca,'XTickLabel',[]);
% set(gca,'XTick',[]);
% set(gca,'YTick',[]);
title('Integrated Spectrum Test','FontSize',20);
ylabel('Cumulative Density','FontSize',18)
xlabel('Averaged Spectrum (1/Hz)','FontSize',18)
% a = get(gca,'XTickLabel');
% set(gca,'XTickLabel',a,'fontsize',16)
axis tight
box off
%%% Plot residual test
subplot(3,2,6)
data  = model_true.data;

nelectrodes = size(data,1);
model_order = model_true.estimated_model_order;
datap = data(:,model_order+1:end);
yestimate = model_spline.signal_estimate;

residuals = datap-yestimate;
autocorr(residuals);
title('Autocorrelation of Residuals','FontSize',20);
xlabel('Lag (s)','FontSize',18);
ylabel('Autocorrelation','FontSize',18);
set(gca,'XTickLabel',[0 model_true.taxis(1:19)],'FontSize',16)
box off
%%
%%%% Fig 2
Labels = {'Standard', 'Spline'};
figure;

%%%%% Plot bar Comp Times
subplot 131
barplot(Labels,ct_standard(2:end),ct_spline(2:end))
ylabel('Computation time (s)','FontSize',20)
 set(gca,'xlim',[0.5 2.5])
% a = get(gca,'YTickLabel');
% set(gca,'YTickLabel',a,'fontsize',16)
 box off
%  xlim =([.75 2.25]);
axis tight

%%% Plot bar GR test stat
subplot 132
barplot(Labels,ts_stand,ts_spline);
xlim=[.5 2.5];
hold on
plot(xlim,[2.2414 2.2414],'--r','LineWidth',2.5)
ylabel('Grenander & Rosenblatt Statistic','FontSize',20)
axis tight
box off
a = get(gca,'YTickLabel');
set(gca,'YTickLabel',a,'fontsize',16)

%%% Plot bar DW stat
subplot 133
barplot(Labels,dwstandard,dwspline)
hold on
plot(xlim,[1 1],'--r','LineWidth',2.5)
plot(xlim,[3 3],'--r','LineWidth',2.5)
ylabel('Durbin-Watson Statistic','FontSize',20)
axis tight
box off
a = get(gca,'YTickLabel');
set(gca,'YTickLabel',a,'fontsize',16)

%%% Save ----------

save('fig1')

% 
 h = get(0,'children');
 for i=1:length(h)
     saveas(h(i), ['fig1_single_node'  num2str(i) 'lowfreq'], 'fig');
     saveas(h(i), ['fig1_single_node'  num2str(i) 'lowfreq'], 'jpg');

 end
 close all;



