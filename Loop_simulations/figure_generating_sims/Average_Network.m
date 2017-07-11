%%%% 
ntrials=20;

config_spline;

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
mnstandard = mean(standardtime);
mnspline = mean(splinetime);

semstandard = std(standardtime)/sqrt(ntrials);
semspline = std(splinetime)/sqrt(ntrials);

figure;
subplot 121
bar([mean(standardtime),mean(splinetime)],'FaceAlpha',0.2)
hold on;
plot([1 1], [mnstandard-2*semstandard, mnstandard+2*semstandard],'k','LineWidth',2);
plot([2 2], [mnspline-2*semspline, mnspline+2*semspline],'k','LineWidth',2);
hold off;
colormap(gray)
Labels = {'Standard', 'Spline'};
set(gca, 'XTick', 1:2, 'XTickLabel', Labels);
ylabel('Computation time (seconds)')



mnstandard = mean(accstand);
mnspline = mean(accspline);

semstandard = std(accstand)/sqrt(ntrials);
semspline = std(accspline)/sqrt(ntrials);
subplot 122
bar([mean(accstand),mean(accspline),],'FaceAlpha',0.2)
hold on;
plot([1 1], [mnstandard-2*semstandard, mnstandard+2*semstandard],'k','LineWidth',2);
plot([2 2], [mnspline-2*semspline, mnspline+2*semspline],'k','LineWidth',2);
hold off;
colormap(gray)
Labels = {'Standard', 'Spline'};
set(gca, 'XTick', 1:2, 'XTickLabel', Labels);
ylabel('Accuracy')
  

figure;
subplot 211
plotchannels(model_true.taxis,model_true.data')
list=prism(9);
subplot 212
for i =1:9
    mySpec(model_true.data(i,:),model_true.sampling_frequency,'yesplot','tapers',list(i,:));
    hold on;
    
end

