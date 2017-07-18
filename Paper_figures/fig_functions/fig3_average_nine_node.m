%%%%
ntrials=100;
config_spline;
model_true.true_coefficients = nine_node_order20_rdi; %%%% MODIFY COEFFICIENTS HERE!
model_true.model_coefficients = model_true.true_coefficients;

nelectrodes = size(model_true.model_coefficients,1);

trials_spline = zeros(nelectrodes,nelectrodes,ntrials);
trials_stand = zeros(nelectrodes,nelectrodes,ntrials);

splinetime = zeros(1,ntrials);
standardtime = zeros(1,ntrials);

for k = 1:ntrials
    
    simulate_network;
    
    tic
    [ adj_spline] = build_ar_splines( model_true);
    splinetime(k)  = toc;
    
    tic
    [ adj_mat] = build_ar( model_true);
    standardtime(k)  = toc;
    
    trials_spline(:,:,k) = adj_spline;
    trials_stand(:,:,k) = adj_mat;
    fprintf([num2str(k) ' \n']);
end

avg_network = mean(trials_spline,3);
std_network = std(trials_spline,0,3);
b = model_true.model_coefficients;

%%% Figure 1
figure;
%%% plot coefficient strengths
subplot 221
[bb, ii] = max(abs(b),[],3);
plotNetwork(bb)
title('Coefficient Strength','FontSize',20)
colorbar
caxis([0 1])

%%% plot one instance of spline inference
subplot 222
plotNetwork(adj_spline)
title('Spline Network','FontSize',20)
colorbar
caxis([0 1])

%%% plot average network
subplot 223
plotNetwork(avg_network);
title('Averaged Network','FontSize',20);
colorbar
caxis([0 1])

%%% plot std network
subplot 224
plotNetwork(std_network)
title('Standard Deviation','FontSize',20);
colorbar
caxis([0 1])

adj_true = model_true.true_coefficients;
adj_true(adj_true~=0)=1;
adj_true=sum(adj_true,3);
adj_true(adj_true~=0)=1;
for i = 1:ntrials
    accspline(i) = network_accuracy(adj_true,trials_spline(:,:,i));
    accstand(i) = network_accuracy(adj_true,trials_stand(:,:,i));
    
end

%%% Figure 2
Labels = {'Standard', 'Spline'};
figure;
%%% computational time bar plot
subplot 121
barplot(Labels, standardtime,splinetime)
set(gca,'xlim',[0.5 2.5])
ylabel('Computation time (s)','FontSize',18)

%%% accuracy bar plot
subplot 122
barplot(Labels,accstand,accspline)
ylabel('Accuracy','FontSize',18)
set(gca,'xlim',[0.5 2.5])
ylim([0 1])


%%% Figure 3
figure;

%%% plot signals
subplot 211
plotchannels(model_true.taxis,model_true.data')

ylabel('Signal','FontSize',18);
xlabel('Time (s)','FontSize',18);
title('Nine Node Network Simulation','FontSize',20);

%%% plot spectrogram
subplot 212
list=lines(9);
for i =1:9
    mySpec(model_true.data(i,:),model_true.sampling_frequency,'yesplot','tapers',list(i,:));
    hold on;
    
end
ylim([-100,0]);
title('Spectrogram','FontSize',20);
ylabel('Power (dB)','FontSize',18)
xlabel('Frequecy (Hz)','FontSize',18)
box off
a = get(gca,'YTickLabel');
set(gca,'YTickLabel',a,'fontsize',16)

h = get(0,'children');
for i=1:length(h)
    
    saveas(h(i), ['fig3_avgntwk'  num2str(i)], 'fig');   
end
close all;


