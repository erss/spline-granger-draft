function [ y] = myPrediction( data, modelCoef, nlags)
% MYPREDICTION computes the one-step prediction for an AR network model of
% coefficients
%
% INPUTS:
%  data           = A matrix of electode data with dimensions electrodes x
%                    time
%  modelCoef      = model coefficients with dimensions electrodes x number 
%                    predictor variables
%  nlags          = The number of lags used as used for each signal
% 
% OUTPUTS:
%  y =  The one-step prediction for each electrode in network with
%       dimensions electrodes x 1.

nelectrodes = size(data,1);

%%% Build predictor variables--------------------------------------------
   X= [];
for k = 1:nelectrodes
    X_temp = []; 
    sgnl = data(k,:)';
    for i=1:nlags                                   %For each lag,
        X_temp = [X_temp, circshift(sgnl,i)];   %... shift x and store it.
     
    end
    X_temp = X_temp(1,:); % use most recent data 
    X = [X X_temp];
end

%%% Put into model--------------------------------------------------------
  y = X*modelCoef';
  y = y';

end

