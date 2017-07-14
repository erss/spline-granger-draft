%%%% 
ntrials=20;
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
subplot 131
[bb, ii] = max(abs(b),[],3);
plotNetwork(bb)
title('Coefficient Strength')
colorbar
caxis([0 1])
subplot 132
plotNetwork(avg_network);
title('Averaged Network');
colorbar
caxis([0 1])
subplot 133
plotNetwork(std_network)
title('Standard Deviation');
colorbar
caxis([0 1])
%%
figure;
strength = reshape(bb,[1 81]);
averages = reshape(avg_network,[1 81]);
stds = reshape(std_network,[1 81]);

plot(strength,averages,'.','MarkerSize',25);
hold on;
plot(strength,stds,'.','MarkerSize',25);
h=legend('Average','Standard Deviation');
set(h,'FontSize',20)
xlabel('Strength of Influence','FontSize',20);
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

%%
Labels = {'Standard', 'Spline'};
figure;
subplot 121
barplot(Labels, standardtime,splinetime)
ylabel('Computation time (s)','FontSize',15)
subplot 122
barplot(Labels,accstand,accspline)
ylabel('Accuracy','FontSize',15)
  

figure;
subplot 211

plotchannels(model_true.taxis,model_true.data')

ylabel('Signal','FontSize',15);
xlabel('Time (s)','FontSize',15);

title('Nine Node Network Simulation','FontSize',15);

subplot 212

list=lines(9);
for i =1:9
    mySpec(model_true.data(i,:),model_true.sampling_frequency,'yesplot','notapers',list(i,:));
    hold on;
  
end
%     h = get(0,'children');
%     for i=1:length(h)
% 
%             saveas(h(i), ['avgntwk'  num2str(i)], 'fig');
%              saveas(h(i), ['avgntwk'  num2str(i)], 'jpg');
%        
%         
%     end
%     close all;
    

