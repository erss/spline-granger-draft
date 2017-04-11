function [ b, yhat ] = estimate_coef( data, adj_mat, model_order, flag,c_pt_times)
% ESTIMATE_COEF builds an AR model using splines, given a known network configuration.
% DATA has dimension number of electrodes by time, ADJ_MAT contains network
% connectivity, model_order is the number of lags for each electrode used to fit
% the model.  B and YHAT are matrices containing the model fit coefficients and signal, 
% respectively, for each electrode in each row.
% flag =0; general AR
% flag = 1; splines


   nelectrodes = size(data,1);            % number electrodes
   nobservations = length(data(1,model_order+1:end)); % number of observations

%b = zeros(nelectrodes,model_order*nelectrodes);
b = zeros(nelectrodes,nelectrodes,model_order);
yhat = zeros(nelectrodes,size(data,2));
yhat = yhat(:,model_order+1:end);
%%% Define control points and build predictors
if nargin == 4
    c_pt_times = [0:10:model_order] ;  % Define Control Point Locations
end

            
s = 0.5;                                    % Define Tension Parameter

% Construct spline regressors for case nelectrodes = 1.
c_pt_times_all = [c_pt_times(1)-2 c_pt_times c_pt_times(end)+2];
Z = zeros(model_order,length(c_pt_times_all));
num_c_pts = length(c_pt_times_all);  %number of control points in total
for i=1:model_order
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
    
    if sum(data_copy) ~= 0
    Z1 = kron(eye(size(data_copy,1)),Z);     % Build block diagonal spline regressors
  
          X = [];                                 % Build history matrix
    for k = 1:size(data_copy,1)
        X_temp = []; 
        sgnl = data_copy(k,:)';
        for i=1:model_order                                   %For each lag,
            X_temp = [X_temp, circshift(sgnl,i)];   %... shift x and store it.
        end
        X_temp = X_temp(model_order+1:end,:);  
        X = [X X_temp];
    end
    
    
    %%% Build Model
     
    % Generate observations for given y
        y = data(ii,model_order+1:end);   
        y = y';

        
        if flag == 1 % splines
    % Fit full model and calculate RSS
       Xfull = X * Z1;      % regressors for y_hat = X*Z1*alpha
       [alpha,~,~] = glmfit(Xfull,y,'normal','constant','off');  % estimate values at control points, alpha
     %  bhat = Z1*alpha(2:end);             % calculate beta values, for every point in space
         bhat = Z1*alpha;                                     % only for electrodes in network
       
         yhat(ii,:) = glmval(alpha,Xfull,'identity','constant','off'); % Get signal estimate.
         
        
        end
        
        if flag == 0  % traditional AR
          [bhat,~,~] = glmfit(X,y,'normal','constant','off');
          yhat(ii,:) = glmval(bhat,X,'identity','constant','off'); % Get signal estimate.

        end
                                           
%         j =1;
%         for k = 1:nelectrodes             
%            if preds(k)
%                b(ii,((k-1)*model_order + 1): k*model_order) = bhat(((j-1)*model_order + 1): j*model_order);         
%                j = j+1;
%            end
%         end
        
       
        
   
    end

end
 b=b;   
%b = reshape(b,[nelectrodes nelectrodes model_order]);
end

