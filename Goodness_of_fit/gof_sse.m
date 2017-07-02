function [ RSQADJ ] = gof_sse( model_spline )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here


data = model_spline.data;
data = data(:,model_spline.estimated_model_order+1:end);
yhat = model_spline.signal_estimate;
b = model_spline.model_coefficients;

residuals = data - yhat;

r = size(b,1)*size(b,3); %CHECK  % number of regressors in model
s = size(residuals,2);
SSR = sum(residuals.^2,2)';          % residuals sum of squares
SST = sum(data.^2,2)';          % total sum of squares

RSQ = 1 - SSR./SST;

RSQADJ = 1 - ((s-1)/(s-r-1))*(1-RSQ);

end

