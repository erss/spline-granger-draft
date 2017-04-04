function [ y] = myPrediction( data, model_coefficients, nlags)
% MYPREDICTION computes the one-step prediction for an AR network model of
% coefficients
%
% INPUTS:
%  data           = A matrix of electode data with dimensions electrodes x
%                    time
%  model_coef      = model coefficients with dimensions electrodes x number 
%                    predictor variables
%  nlags          = The number of lags used as used for each signal
% 
% OUTPUTS:
%  y =  The one-step prediction for each electrode in network with
%       dimensions [number of electrodes] x [1].
%
% Example for simulating a signal:
% for k = nlags:length(data)-1;
%     data(:,k+1) = myPrediction(data(:,1:k),model_coef,nlags);
%     data(:,k+1) = data(:,k+1) + noise.*randn(size(data,1),1);
% end
%
% NOTE model_coefficents and nlags depend on each other

nelectrodes = size(data,1);

%%% Build predictor variables--------------------------------------------
   
X= [];
   
for n = 1:nelectrodes
    X = [X, data(n,end:-1:end-nlags+1)];
end

%%% Put into model--------------------------------------------------------
  y = X*model_coefficients';
  y = y';

end

