%%%%
ntrials=10;
config_spline;
model_true.true_coefficients = nine_node_order20_rdi; %%%% MODIFY COEFFICIENTS HERE!
model_true.model_coefficients = model_true.true_coefficients;


%%
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

figure;
subplot 221
[bb, ii] = max(abs(b),[],3);
plotNetwork(bb)
title('Coefficient Strength','FontSize',20)
colorbar
caxis([0 1])

subplot 222
plotNetwork(adj_spline)
title('Spline Network','FontSize',20)
colorbar
caxis([0 1])

subplot 223
plotNetwork(avg_network);
title('Averaged Network','FontSize',20);
colorbar
caxis([0 1])
subplot 224
plotNetwork(std_network)
title('Standard Deviation','FontSize',20);
colorbar
caxis([0 1])
%%
% figure;
% strength = reshape(bb,[1 81]);
% averages = reshape(avg_network,[1 81]);
% stds = reshape(std_network,[1 81]);
%
% plot(strength,averages,'.','MarkerSize',25);
% hold on;
% plot(strength,stds,'.','MarkerSize',25);
% h=legend('Average','Standard Deviation');
% set(h,'FontSize',20)
% xlabel('Strength of Influence','FontSize',20);
%%
% bb= sum(b,3);
% bb(bb~=0) = 1;
% figure;
% plotNetwork(bb)



%
adj_true = model_true.true_coefficients;
adj_true(adj_true~=0)=1;
adj_true=sum(adj_true,3);
adj_true(adj_true~=0)=1;
for i = 1:ntrials
    accspline(i) = network_accuracy(adj_true,trials_spline(:,:,i));
    accstand(i) = network_accuracy(adj_true,trials_stand(:,:,i));
    
    
end

%%% Figure 1
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

figure;
subplot 211

plotchannels(model_true.taxis,model_true.data')
% hold on;
%  plot([0,.2], [min(model_true.data(1,:)), min(model_true.data(1,:))], 'r', 'LineWidth', 2.5);


ylabel('Signal','FontSize',18);
xlabel('Time (s)','FontSize',18);
title('Nine Node Network Simulation','FontSize',20);

subplot 212
list=lines(9);
for i =1:9
    mySpec(model_true.data(i,:),model_true.sampling_frequency,'yesplot','tapers',list(i,:));
    hold on;
    
end
ylim([-120,0]);
box off
%     h = get(0,'children');
%     for i=1:length(h)
%
%             saveas(h(i), ['avgntwk'  num2str(i)], 'fig');
%              saveas(h(i), ['avgntwk'  num2str(i)], 'jpg');
%
%
%     end
%     close all;


