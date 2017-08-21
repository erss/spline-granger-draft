%%%%% single node analysis
%%% set model coefficitens to single_node_order20 & single_node_low_freq
%clear all;
ntrials = 1000;
ct_spline = zeros(1,ntrials);
ct_standard = zeros(1,ntrials);
ts_spline = zeros(1,ntrials);
ts_stand = zeros(1,ntrials);
dwstandard = zeros(1,ntrials);
dwspline = zeros(1,ntrials);
fails =[];
fails_st =[];

for i = 1:ntrials
    config_spline;
    model_true.noise_type = 'white';
    model_true.true_coefficients = single_node_order20; %%%% MODIFY COEFFICIENTS HERE!
    
    
    model_true.model_coefficients = model_true.true_coefficients;
    simulate_network;
    infer_network;
    
    ct_spline(i)  = model_spline.computation_time;
    ct_standard(i) = model_standard.computation_time;
    
    
    acc_spline(i)  = model_spline.accuracy;
    acc_standard(i) = model_standard.accuracy;
    
    [nw, dwstandard(i)] = dwstat(model_standard);
    fails_st = [fails nw];
    [nw, dwspline(i)] = dwstat(model_spline);
    fails = [fails nw];
    
    
      [m2fit, m3fit] = grstat1(model_true,model_spline,model_standard);
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

%%% Plot spectrogram
subplot(2,2,2)

mySpec(model_true.data(1,:),model_true.sampling_frequency,'yesplot','tapers');
box off
axis tight
a = get(gca,'YTickLabel');
set(gca,'YTickLabel',a,'fontsize',16)
title('Spectrogram','FontSize',20);
ylabel('Power (dB)','FontSize',18)
xlabel('Frequency (Hz)','FontSize',18)
%%% Plot coefficients
subplot(2,2,[3 4])
gof_bootstrap(model_true,model_spline,model_standard);
box off
title('AR Coefficients','FontSize',20);
ylabel('Magnitude','FontSize',18)
xlabel('Lag (s)','FontSize',18)
a = get(gca,'XTickLabel');
set(gca,'XTickLabel',a,'fontsize',16)
plot([0 0.06],[0 0],'color',[.57 .57 .57],'LineWidth',1.7)

    h = legend('Spline Estimated','Standard Estimated','True');
        set(h,'FontSize',15,'Location','NorthEast');% 
%%%% Fig 2
Labels = {'Standard', 'Spline'};
figure;

%%%%% Plot bar Comp Times
subplot 131
barplot(Labels,ct_standard,ct_spline)
ylabel('Computation time (s)','FontSize',20)
set(gca,'xlim',[0.5 2.5])
box off
a = get(gca,'YTickLabel');
set(gca,'YTickLabel',a,'fontsize',16)

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
% %%% Save
% h = get(0,'children');
% for i=1:length(h)
% saveas(h(i), ['fig2_single_node'  num2str(i) 'specpeak'], 'fig');     
% end
% close all;
% 
% 
% 
