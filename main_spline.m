clear all;
close all;
clc;
%% Configure parameters

config_spline;

%%% Simulate data
simulate_network;
% mySpec(model_true.data(1,:),model_true.sampling_frequency,'yesplot','tapers'); 
%%
%%% Fit Networks
infer_network;

%%%% Plot stuff
plot_data;
%%%%

%%
% gof_bootstrap(model_true,model_spline,model_standard);
% gof_residuals(model_spline)
% gof_spectrum(model_true,model_spline,model_standard);
% figure; plotchannels(model_true.data');
 %for i =1:5; figure; mySpec(model_true.data(i,:),f0,'yesplot','tapers'); end

%% Goodness of Fit & Plots
% order = model_standard.estimated_model_order;
% data = model_standard.data;
% splinedata =   model_standard.signal_estimate; 
% [ cons1 ] = gof_consistency( model_standard)  %%% compute consistency of inferred network
% consistency(data,data(:,order+1:end)-splinedata)

% gof_residuals(model_spline)
%   model_spline.model_coefficients=single_node_low_freq; %%% change coefficients
% model_spline.estimated_model_order = 2; %%% computre consistency of 'incorrect' network
%  [ cons2] = gof_consistency( model_spline);



%[ sse ] = gof_sse( model_spline)

% gof_derbin_h( model_true, model_spline)
 gof_bootstrap(model_true,model_spline,model_standard);
% gof_residuals(model_spline);
% gof_spectrum(model_true,model_spline,model_standard);

% %% Compute g.o.f with Bartlett test.
% Ij_low = [0,4; 4,8;  8,12; 12,20; 20,30; 30,50];
% [result] = gof_bartlett(model_spline, Ij_low);
%       plot_gof_bartlett(result,       Ij_low);
% 
% Ij_high = (50:5:120)';  Ij_high = [Ij_high(1:end-1), Ij_high(2:end)];
% [result] = gof_bartlett(model_spline, Ij_high);
%       plot_gof_bartlett(result,       Ij_high);
