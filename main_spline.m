clear all;
close all;
clc;
%% Configure parameters

config_spline;

%%% Simulate data
simulate_network;
%%
%%% Fit Networks
infer_network;

%%%% Plot stuff
plot_data;
%%%%

%% Goodness of Fit & Plots
goodness_of_fit;

