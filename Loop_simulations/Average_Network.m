%%%% 
n=20;

config_spline;


nelectrodes = size(model_true.model_coefficients,1);

trials = zeros(nelectrodes,nelectrodes,n);
trials_stand = zeros(nelectrodes,nelectrodes,n);

splinetime = zeros(1,n);
standardtime = zeros(1,n);

for k = 1:n
    
   simulate_network;
   
   tic
[ adj_spline] = build_ar_splines( model_true);
splinetime(k)  = toc;

   tic
[ adj_mat] = build_ar( model_true);
standardtime(k)  = toc;
   
   trials(:,:,k) = adj_spline;
 trials_stand(:,:,k) = adj_mat;
    fprintf([num2str(k) ' \n']);
end

avg_network = mean(trials,3);
std_network = std(trials,0,3);

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
for i = 1:n
    distance_spline(i) = jdist(adj_true,trials(:,:,i));
    distance_stand(i) = jdist(adj_true,trials_stand(:,:,i));
    
    
end

figure;
subplot 121
bar([mean(standardtime),mean(splinetime)])
colormap(gray)
Labels = {'Standard', 'Spline'};
set(gca, 'XTick', 1:2, 'XTickLabel', Labels);
ylabel('Computation time (seconds)')

subplot 122
bar([1-mean(distance_stand),1-mean(distance_spline),])
colormap(gray)
Labels = {'Standard', 'Spline'};
set(gca, 'XTick', 1:2, 'XTickLabel', Labels);
ylabel('Jaccard Similarity')
   
