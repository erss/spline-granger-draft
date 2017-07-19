%%%%% single node analysis
%%% set model coefficitens to single_node_order20 & single_node_low_freq
clear all;
ntrials = 2;
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
    model_true.true_coefficients = single_node_low_freq; %%%% MODIFY COEFFICIENTS HERE!
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
    
    [ts_spline(i), ts_stand(i)] = grstat(model_true,model_spline,model_standard);
 
end


%%% Fig 1
figure;

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
a = get(gca,'XTickLabel');
set(gca,'XTickLabel',a,'fontsize',16)

%%% Plot spectral test
subplot(3,2,5)
gof_spectrum(model_true,model_spline,model_standard);
set(gca,'YTickLabel',[]);
set(gca,'XTickLabel',[]);
set(gca,'XTick',[]);
set(gca,'YTick',[]);
title('Integrated Spectrum Test','FontSize',20);
ylabel('Cumulative Density','FontSize',18)
xlabel('Averaged Spectrum (1/Hz)','FontSize',18)
a = get(gca,'XTickLabel');
set(gca,'XTickLabel',a,'fontsize',16)

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

%%%% Fig 2
Labels = {'Standard', 'Spline'};
figure;

%%%%% Plot bar Comp Times
subplot 131
barplot(Labels,ct_standard,ct_spline)
ylabel('Computation time (s)','FontSize',20)
set(gca,'xlim',[0.5 2.5])
a = get(gca,'YTickLabel');
set(gca,'YTickLabel',a,'fontsize',16)
box off
%  xlim =([.75 2.25]);

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

%%% Save
% h = get(0,'children');
% for i=1:length(h)
%     saveas(h(i), ['fig1_single_node'  num2str(i) 'lowfreq'], 'fig');
%     
%     
%     
% end
% close all;



