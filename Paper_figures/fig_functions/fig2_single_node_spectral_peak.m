%%%%% single node analysis
%%% set model coefficitens to single_node_order20 & single_node_low_freq

clear all;
ntrials = 1; %1000
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
    model_true.true_coefficients =single_node_order20;% %%%% MODIFY COEFFICIENTS HERE!
    model_true.model_coefficients = model_true.true_coefficients;   

%%% Define model inputs for spline Granger & standard Granger -------------

model_true.s = 0.5;                     % tension parameter for spline
model_true.estimated_model_order = 30;  % model_order used to estimate

model_true.cntrl_pts = [0:5:model_true.estimated_model_order];
%%% Define network testing parameters -------------------------------------

model_true.q = 0.05;            % FDR max number acceptable proportion of false discoveries
model_true.nsurrogates = 1000;   % number of surrogates used for bootstrapping
model_true.nrealizations = 20; % number of realizations used for spectral testing



for i = 1:ntrials
%%%% MODIFY COEFFICIENTS HERE!
    
    simulate_network;
    infer_network;
    
    ct_spline(i)  = model_spline.computation_time;
    ct_standard(i) = model_standard.computation_time;
    
    
    acc_spline(i)  = model_spline.accuracy;
    acc_standard(i) = model_standard.accuracy;
    
    [nw, dwstandard(i)] = dwstat(model_standard);
    fails_standard = [fails_standard nw];
    [nw, dwspline(i)] = dwstat(model_spline);
    fails = [fails nw];
    
    
[m2fit] = grstat(model_true,model_spline)
    [m3fit] = grstat(model_true,model_standard)
 ts_spline(i) = m2fit.stat;
 ts_stand(i)  = m3fit.stat;

end

%% Fig 1
figure;

%%% Plot signal trace
subplot(2,2,1)
plot(model_true.taxis,model_true.data(1,:),'k', 'LineWidth', 1);
hold on;   
plot([0,.2], [min(model_true.data(1,:)), min(model_true.data(1,:))], 'k', 'LineWidth', 2.5);
box off
title('Simulated Signal','FontSize',20);
xlabel('Time (s)','FontSize',18);
a = get(gca,'YTickLabel');
set(gca,'YTickLabel',a,'fontsize',16)
box off
set(gca,'YTickLabel',[]);
set(gca,'XTickLabel',[]);
set(gca,'XTick',[]);
set(gca,'YTick',[]);
%%
%%% Plot spectrogram
subplot(2,2,2)

[faxis, Sxx] =mySpec(model_true.data(1,:),model_true.sampling_frequency,'tapers');
        plot((faxis),10*log(Sxx),'col','k','LineWidth',2);
xlim=([0 model_true.sampling_frequency/4]);

title('Spectrogram','FontSize',20);
ylabel('Power (dB)','FontSize',18)
xlabel('Frequency (Hz)','FontSize',18)
box off
axis tight


%%% Plot coefficients
subplot(2,2,[3 4])
gof_bootstrap(model_true,model_spline,model_standard);
box off
title('Estimated Model Coefficients','FontSize',20);
ylabel('Magnitude','FontSize',18)
xlabel('Lag (s)','FontSize',18)
hold on; 
plot([0 0.06],[0 0],'color',[.57 .57 .57],'LineWidth',1.7)

    h = legend('Spline Estimated','Standard Estimated','True');
        set(h,'FontSize',15,'Location','NorthEast');% 
%%%% Fig 2
Labels = {'Standard', 'Spline'};
figure;

%%%%% Plot bar Comp Times
subplot 131
barplot(Labels,ct_standard(2:end),ct_spline(2:end))
xlim=[.5 2.5];
hold on
plot(xlim,[0 0])
ylabel('Computation time (s)','FontSize',20)

box off
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


%%% Plot bar DW stat
subplot 133
barplot(Labels,dwstandard,dwspline)
hold on
plot(xlim,[1 1],'--r','LineWidth',2.5)
plot(xlim,[3 3],'--r','LineWidth',2.5)
ylabel('Durbin-Watson Statistic','FontSize',20)
axis tight
box off


%% %%% Save
% 
% save('fig2')
% 
% h = get(0,'children');
% for i=1:length(h)
% saveas(h(i), ['fig2_single_node'  num2str(i) 'specpeak'], 'fig');     
% saveas(h(i), ['fig2_single_node'  num2str(i) 'specpeak'], 'jpg');     
% 
% end
% close all;



