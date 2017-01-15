function [ b ] = estimate_coef( data, adj_mat, nlags)
% ESTIMATE_COEF builds an AR model using splines, given a known network configuration.
% DATA has dimension number of electrodes by time, ADJ_MAT contains network
% connectivity, NLAGS is the number of lags for each electrode used to fit
% the model.  B is a matrix containing the model fit coefficients for each
% electrode in each row.

   nelectrodes = size(data,1);            % number electrodes
   nobservations = length(data(1,nlags+1:end)); % number of observations

b = zeros(nelectrodes,nlags*nelectrodes);
%%% Define control points and build predictors

c_pt_times = [0:10:nlags];  % Define Control Point Locations

            
s = 0.5;                                    % Define Tension Parameter

% Construct spline regressors for case nelectrodes = 1.
c_pt_times_all = [c_pt_times(1)-2 c_pt_times c_pt_times(end)+2];
Z = zeros(nlags,length(c_pt_times_all));
num_c_pts = length(c_pt_times_all);  %number of control points in total
for i=1:nlags
    nearest_c_pt_index = max(find(c_pt_times_all<i));
    nearest_c_pt_time = c_pt_times_all(nearest_c_pt_index);
    next_c_pt_time = c_pt_times_all(nearest_c_pt_index+1);
    u = (i-nearest_c_pt_time)/(next_c_pt_time-nearest_c_pt_time);
    p=[u^3 u^2 u 1]*[-s 2-s s-2 s;2*s s-3 3-2*s -s;-s 0 s 0;0 1 0 0];
    Z(i,:) = [zeros(1,nearest_c_pt_index-2) p zeros(1,num_c_pts-4-(nearest_c_pt_index-2))];
end



%%% Build model for every electrode in network

for ii = 1:nelectrodes
    % Build history regressors
   
    preds = logical(adj_mat(ii,:)); %% Use results of F-test to input in network
    data_copy = data(preds,:); %% Remove electrodes not connected in spline network
    
    
    Z1 = kron(eye(size(data_copy,1)),Z);     % Build block diagonal spline regressors
  
     X = [];                                 % Build history matrix
    for k = 1:size(data_copy,1)
        X_temp = []; 
        sgnl = data_copy(k,:)';
        for i=1:nlags                                   %For each lag,
            X_temp = [X_temp, circshift(sgnl,i)];   %... shift x and store it.
        end
        X_temp = X_temp(nlags+1:end,:);  
        X = [X X_temp];
    end
    
    
    %%% Build Model
     
    % Generate observations for given y
        y = data(ii,nlags+1:end);   
        y = y';

    % Fit full model and calculate RSS
       Xfull = X * Z1;      % regressors for y_hat = X*Z1*alpha
       [alpha,~,stats] = glmfit(Xfull,y);  % estimate values at control points, alpha
       bhat = Z1*alpha(2:end);             % calculate beta values, for every point in space
                                           % only for electrodes in network
         
        j =1;
        for k = 1:nelectrodes             
           if preds(k)
               b(ii,((k-1)*nlags + 1): k*nlags) = bhat(((j-1)*nlags + 1): j*nlags);         
               j = j+1;
           end
        end

    end

end

