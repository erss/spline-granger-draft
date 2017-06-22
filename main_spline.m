clear all; close all;
clc;
%% Configure parameters
config_spline;

%% Simulate data
simulate_network;

%% Fit Networks
infer_network;

%% Plot stuff
plot_data;
% figure; plotchannels(model_true.data'); for i =1:3; figure; mySpec(model_true.data(i,:),f0,'yesplot','tapers'); end

%% Goodness of Fit & Plots

% gof_bootstrap(model_true,model_spline,model_standard);
% gof_residuals(model_spline);
% gof_spectrum(model_true,model_spline,model_standard);