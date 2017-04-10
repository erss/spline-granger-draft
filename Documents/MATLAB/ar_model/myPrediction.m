function [ y] = myPrediction( data, model_coefficients)
% MYPREDICTION computes the one-step prediction for an AR network model of
% coefficients
%
% INPUTS:
%  data           = A matrix of electode data with dimensions electrodes x
%                    time
%  model_coef      = model coefficients with dimensions electrodes x
%                              electrodes x number time lags, i.e. if
%                              model_coeff is 3 x 3 x 10, then there are 3
%                              electrodes in the network whose activity is
%                              deterined by 10 time steps in the past.  The
%                              influence at lag = 1 is given by 
%                              model_coeff(:,:,1). This gives the network
%                              connectivity at a lag 1, where the rows
%                              represent each electrode fit, i. e.
%                              model_coeff(1,3,1) is the influence of
%                              electrode 3 on electrode 1 at one time lag
%                              in the past
%                             
% OUTPUTS:
%  y =  The one-step prediction for each electrode in network with
%       dimensions [number of electrodes] x [1].
%
% Example for simulating a signal:
% for k = nlags:length(data)-1;
%     data(:,k+1) = myPrediction(data(:,1:k),model_coef);
%     data(:,k+1) = data(:,k+1) + noise.*randn(size(data,1),1);
% end
%

nelectrodes = size(data,1);
nlags = size(model_coefficients,3);

%%% Build predictor variables--------------------------------------------
  
for p = 1:nelectrodes
    X = [];
    b = [];  
    for n = 1:nelectrodes
        X = [X, data(n,end:-1:end-nlags+1)];
        b = [b, squeeze(model_coefficients(p,n,:))'];

    end
     y(p) = X*b';
  
end
%%% Put into model--------------------------------------------------------

  y = y';

end

