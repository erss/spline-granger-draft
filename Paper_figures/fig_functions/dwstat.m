function [notwhite, dw,pval] = dwstat( model)
%UNTITLED8 Summary of this function goes here
%   Detailed explanation goes here
%%%% goodness_of_fit_residuals
data  = model.data;
nelectrodes = size(data,1);
model_order = model.estimated_model_order;
datap = data(:,model_order+1:end);
yestimate = model.signal_estimate;
dt = 1/model.sampling_frequency;
T = model.T;
taxis = dt:dt:T';


residuals = zeros(nelectrodes,length(datap));

for electrode = 1:nelectrodes
    true_signal =  datap(electrode,:);
    est_signal = yestimate(electrode,:);

    residuals(electrode,:) = true_signal-est_signal;
end

 [dw, pval] = whiteness(data,residuals);
% % A standard rule of thumb is that |dw < 1| or |dw > 3| indicates a high
% % chance of residuals serial correlation; this implies poor VAR model fit.
 sig = significance(pval,0.05,'FDR')
notwhite = find(sig);
if isempty(notwhite)
    fprintf('all residuals are white by Durbin-Watson test at significance %g\n',0.05);
else
    fprintf(2,'WARNING: autocorrelated residuals at significance %g for variable(s): %s\n',0.05,num2str(notwhite));
end


end

