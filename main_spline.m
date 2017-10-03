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
goodness_of_fit;

