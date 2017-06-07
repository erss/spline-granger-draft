function [bhat, yhat] = estimate_standard( data, adj_mat, model_order )
% ESTIMATE_COEFFICIENT_FITS builds network model from MVAR spline modeling 
%
% INPUTS:
%  data           = A matrix of electode data with dimensions electrodes x
%                    time
%  adj_mat        = network configuration
%  model_order         = The number of lags used as used for predictor variables
%  c_pt_times     = The number of lags used as used for predictor variables
% 
% OUTPUTS:
%  bhat   = coefficient estimates
%  yhat    = signal estimates


%% Initialize variables & outputs
    nelectrodes = size(data,1);            % number electrodes
    nobservations = length(data(1,model_order+1:end)); % number of observations
    bhat = zeros(nelectrodes,nelectrodes,model_order);
    yhat = zeros(nelectrodes,size(data,2));
    yhat = yhat(:,model_order+1:end);
%% Define control points and build predictors

for electrode = 1:nelectrodes
    
   
    model_size = sum(adj_mat(electrode,:));
    if model_size > 0
         % Nested model spline regressors

        % Generate observations for given y
        x = data(electrode,:);
        y = x(model_order+1:end);   
        y = y';

        % Build history matrix
         X = []; 
        % Generate regressors for network structure
         indices = logical(adj_mat(electrode,:));
         data_subset = data(indices,:);

        for k = 1:model_size
            X_temp = []; 
            sgnl = data_subset(k,:)';


            for i=1:model_order                                   %For each lag,
                X_temp = [X_temp, circshift(sgnl,i)];   %... shift x and store it.
            end
            X_temp = X_temp(model_order+1:end,:);  


        X = [X X_temp];

        end
       
       % Xfull = X * Z0;  % regressors for y_hat = X0*Z0*a0
        [b,~,~] = glmfit(X,y,'normal','constant','off');

        yhat(electrode,:) = glmval(b,X,'identity','constant','off'); % Get signal estimate.
       

        j=1;
          for p = 1:nelectrodes
            if adj_mat(electrode,p) == 1
                bhat(electrode,p,:) = b(j:j+model_order-1);
                j= j+model_order;
            end    
         end
    end

end

    
  
end

