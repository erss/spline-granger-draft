clear all;
close all;
clc;
%% Configure parameters

config_spline;

%%% Simulate data
simulate_network;
 %mySpec(model_true.data(1,:),model_true.sampling_frequency,'yesplot','tapers'); 
%%
%%% Fit Networks
infer_network;

%%%% Plot stuff
plot_data;
%%%%

%% Goodness of Fit & Plots
figure;
 gof_bootstrap(model_true,model_spline,model_standard);
 gof_residuals(model_spline);
 gof_spectrum(model_true,model_spline,model_standard);
 
 figure;
  [faxis, Sxx] = mySpec(model_true.data,f0,'noplot','notapers');
  plot(faxis,Sxx)
hold on
  [faxis, Sxx] = mySpec(model_spline.signal_estimate,f0,'noplot','notapers');
  plot(faxis,Sxx,'r')

