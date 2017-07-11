%%%%% single node analysis

ntrials = 10;
fails =[];
for i = 1:ntrials 
config_spline;
simulate_network;
infer_network;

ct_spline(i)  = model_spline.computation_time;
ct_standard(i) = model_standard.computation_time;


acc_spline(i)  = model_spline.accuracy;
acc_standard(i) = model_standard.accuracy;

nw = gof_residuals(model_spline);
fails = [fails nw];

[ts_spline(i), ts_stand(i)] = gof_spectrum(model_true,model_spline,model_standard);
close all;
end
%%
figure;
subplot(3,2,1)
plot(model_true.data(1,:));
subplot(3,2,2)
mySpec(model_true.data(1,:),f0,'yesplot','tapers');
subplot(3,2,[3 4])
 gof_bootstrap(model_true,model_spline,model_standard);
subplot(3,2,[5 6])
gof_spectrum(model_true,model_spline,model_standard);

%%
figure;
subplot 121
bar([mean(ct_standard),mean(ct_spline)])
colormap(gray)
Labels = {'Standard', 'Spline'};
set(gca, 'XTick', 1:2, 'XTickLabel', Labels);
ylabel('Computation time (seconds)')

subplot 122
bar([mean(ts_stand),mean(ts_spline),])
colormap(gray)
Labels = {'Standard', 'Spline'};
set(gca, 'XTick', 1:2, 'XTickLabel', Labels);
ylabel('GoF metric')
   


