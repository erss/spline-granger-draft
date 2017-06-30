clear; close all;
clc;
%%% Configure parameters
config_spline;

%%% Simulate data
simulate_network;

%%% Fit Networks
infer_network;

%%%% Plot stuff
plot_data;

% figure; plotchannels(model_true.data');
% for i =1:3; figure; mySpec(model_true.data(i,:),f0,'yesplot','tapers'); end

%% Goodness of Fit & Plots
% gof_derbin_h( model_true, model_spline)
% gof_bootstrap(model_true,model_spline,model_standard);
% gof_residuals(model_spline);
% gof_spectrum(model_true,model_spline,model_standard);

%%% Compute g.o.f with Bartlett test.
Ij_low = [0,4; 4,8;  8,12; 12,20; 20,30; 30,50];
[result] = gof_bartlett(model_spline, Ij_low);
      plot_gof_bartlett(result,       Ij_low);

Ij_high = (50:5:120)';  Ij_high = [Ij_high(1:end-1), Ij_high(2:end)];
[result] = gof_bartlett(model_spline, Ij_high);
      plot_gof_bartlett(result,       Ij_high);
