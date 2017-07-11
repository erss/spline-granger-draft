%%%%% single node analysis
%%% set model coefficitens to single_node_order20 & single_node_low_freq
ntrials = 20;
fails =[];
for i = 1:ntrials 
config_spline;
simulate_network;
infer_network;

ct_spline(i)  = model_spline.computation_time;
ct_standard(i) = model_standard.computation_time;


acc_spline(i)  = model_spline.accuracy;
acc_standard(i) = model_standard.accuracy;
[nw dwstandard(i)] = gof_residuals(model_standard);
fails = [fails nw];
[nw dwspline(i)] = gof_residuals(model_spline);
fails = [fails nw];

[ts_spline(i), ts_stand(i)] = gof_spectrum(model_true,model_spline,model_standard);
close all;

fprintf(num2str(i));
end
%%
figure;
subplot(3,2,1)
plot(model_true.taxis,model_true.data(1,:),'k');
title('Simulated Trace');
xlabel('Time (s)','FontSize',14);
subplot(3,2,2)
mySpec(model_true.data(1,:),model_true.sampling_frequency,'yesplot','tapers');
subplot(3,2,[3 4])
 gof_bootstrap(model_true,model_spline,model_standard);
subplot(3,2,5)
gof_spectrum(model_true,model_spline,model_standard);
subplot(3,2,6)
data  = model_true.data;

nelectrodes = size(data,1);
model_order = model_true.estimated_model_order;
datap = data(:,model_order+1:end);
yestimate = model_spline.signal_estimate;

residuals = datap-yestimate;
autocorr(residuals);

%%

mnstandard = mean(ct_standard);
mnspline = mean(ct_spline);

semstandard = std(ct_standard)/sqrt(ntrials);
semspline = std(ct_spline)/sqrt(ntrials);


figure;
subplot 131
bar([mean(ct_standard),mean(ct_spline)],'FaceAlpha',0.2)
hold on;
plot([1 1], [mnstandard-2*semstandard, mnstandard+2*semstandard],'k','LineWidth',2);
plot([2 2], [mnspline-2*semspline, mnspline+2*semspline],'k','LineWidth',2);
hold off;
colormap(gray)
Labels = {'Standard', 'Spline'};
set(gca, 'XTick', 1:2, 'XTickLabel', Labels);
ylabel('Computation time (seconds)')


mnstandard = mean(ts_stand);
mnspline = mean(ts_spline);

semstandard = std(ts_stand)/sqrt(ntrials);
semspline = std(ts_spline)/sqrt(ntrials);

subplot 132
bar([mean(ts_stand),mean(ts_spline),],'FaceAlpha',0.2)
hold on;
plot([1 1], [mnstandard-2*semstandard, mnstandard+2*semstandard],'k','LineWidth',2);
plot([2 2], [mnspline-2*semspline, mnspline+2*semspline],'k','LineWidth',2);
hold off

colormap(gray)
Labels = {'Standard', 'Spline'};
set(gca, 'XTick', 1:2, 'XTickLabel', Labels);
xlim=get(gca,'xlim');
hold on
plot(xlim,[2.2414 2.2414],'--r','LineWidth',2)
ylabel('Spectral GoF metric')


mnstandard = mean(dwstandard);
mnspline = mean(dwspline);

semstandard = std(dwstandard)/sqrt(ntrials);
semspline = std(dwspline)/sqrt(ntrials);


subplot 133
bar([mnstandard,mnspline,],'FaceAlpha',0.2)
hold on;
plot([1 1], [mnstandard-2*semstandard, mnstandard+2*semstandard],'k','LineWidth',2);
plot([2 2], [mnspline-2*semspline, mnspline+2*semspline],'k','LineWidth',2);
hold off

colormap(gray)
Labels = {'Standard', 'Spline'};
set(gca, 'XTick', 1:2, 'XTickLabel', Labels);
xlim=get(gca,'xlim');
hold on
plot(xlim,[1 1],'--r','LineWidth',2)
plot(xlim,[3 3],'--r','LineWidth',2)
ylabel('DW GoF metric')
   
   


