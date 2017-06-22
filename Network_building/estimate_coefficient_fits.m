function [bhat, yhat] = estimate_coefficient_fits( model, adj_mat)
% ESTIMATE_COEFFICIENT_FITS builds network model from MVAR spline modeling 
%
% INPUTS:
%  data           = A matrix of electode data with dimensions electrodes x
%                    time
%  adj_mat        = network configuration
%  model_order     = The number of lags used as used for predictor variables
%  c_pt_times     = The number of lags used as used for predictor variables
% 
% OUTPUTS:
%  bhat   = coefficient estimates
%  yhat    = signal estimates
s = model.s;
data= model.data;
model_order = model.estimated_model_order;
cntrl_pts  = model.cntrl_pts;

%% Initialize variables & outputs
    nelectrodes = size(data,1);            % number electrodes
    nobservations = length(data(1,model_order+1:end)); % number of observations
    bhat = zeros(nelectrodes,nelectrodes,model_order);
    yhat = zeros(nelectrodes,size(data,2));
    yhat = yhat(:,model_order+1:end);
%% Define control points and build predictors

% if nargin == 3
%     c_pt_times = unique([0:10:model_order model_order]) ;  % Define Control Point Locations
% end
            
% Construct spline regressors.
c_pt_times_all = [cntrl_pts(1)-2 cntrl_pts cntrl_pts(end)+2];
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



for electrode = 1:nelectrodes
    
   
    model_size = sum(adj_mat(electrode,:));
    if model_size > 0
        Z0 = kron(eye(model_size),Z);   % Nested model spline regressors

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
       % 
        Xfull = X * Z0;  % regressors for y_hat = X0*Z0*a0
        [alpha,~,~] = glmfit(Xfull,y,'normal','constant','off');

        yhat(electrode,:) = glmval(alpha,Xfull,'identity','constant','off'); % Get signal estimate.
        b = Z0*alpha;  

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

