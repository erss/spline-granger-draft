function [ output_args ] = gof_derbin_h( model1 , model2 )
%UNTITLED9 Summary of this function goes here
%   Detailed explanation goes here
%%% Compare spectrum of true signal and of estimated signal
%%%%%% NOTE highly dependent on noise used to generate model
% Define model inputs

model_order = model1.estimated_model_order;
data = model1.data;
spline_data = model2.signal_estimate;
residuals = data(:,model_order+1:end) - spline_data;

nelectrodes = size(model1.data,1); % number of electrodes
design_matrix = model2.design_matrix;
covb = model2.covb;
T = length(residuals(1,:)); % length of residuals
for ii = 1:nelectrodes
    cb= covb(ii);
    [~, dw ] = dwtest(residuals(ii,:),design_matrix)
    
    dh = (1-0.5*dw)*sqrt(T/(1-T*cb))
    dhsig = 1-normcdf(dh,0,1)    %[dh, dhsig] = durbinh (dw, covb(ii), length(residuals),2)
end

end

