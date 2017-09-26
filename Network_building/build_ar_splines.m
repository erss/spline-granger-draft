function [ adj_mat] = build_ar_splines( model )
% BUILD_AR_SPLINES builds network model from MVAR modeling and uses
% regression splines to reduce dimensionality.
%
% INPUTS:
%  data           = A matrix of electode data with dimensions electrodes x
%                    time
%  model_order    = The number of lags used as used for predictor variables
%  cntrl_pts      = (OPTIONAL) The number of lags used as used for predictor variables
% 
% OUTPUTS:
%  adj_mat = adjacencey matrix for corresponding network

warning off

data = model.data;
model_order = model.estimated_model_order;
cntrl_pts = model.cntrl_pts;
q = model.q;
s = model.s;
%% Initialize variables & outputs
    nelectrodes = size(data,1);            % number electrodes
    adj_mat = zeros(nelectrodes);
    F = zeros(nelectrodes);                % matrix of F statistics generated by removing subset of
                                           %    predictor variables
    nobservations = length(data(1,model_order+1:end)); % number of observations
%% Define control points and build predictors
            
% Construct spline regressors.
c_pt_times_all = [cntrl_pts(1)-100 cntrl_pts cntrl_pts(end)+100];
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
Z0 = kron(eye(nelectrodes-1),Z);   % Nested model spline regressors
Z1 = kron(eye(nelectrodes),Z);     % Full model spline regressors

% Build history matrix
X = []; 

% Vector of electrode names corresponding to order in build matrix
% e.g. for trivariate case with two lags, [1 1 2 2 3 3]
e_names = zeros(1,model_order*nelectrodes); 
j = 1;
for k = 1:nelectrodes
    X_temp = []; 
    sgnl = data(k,:)';
    for i=1:model_order                                   %For each lag,
        X_temp = [X_temp, circshift(sgnl,i)];   %... shift x and store it.
        e_names(j) = k;
        j = j+ 1;
    end
    X_temp = X_temp(model_order+1:end,:);  
    X = [X X_temp];
end

%% Determine connectivity using F test

% Build models and test correlation for every electrode pair
    for electrode = 1:nelectrodes
        
        % Generate observations for given y
        x = data(electrode,:);
        y = x(model_order+1:end);   
        y = y';


        % Fit full model and calculate RSS
        Xfull = X * Z1;      % regressors for y_hat = X*Z1*alpha
        [alpha,~,stats] = glmfit(Xfull,y,'normal','constant','off');
        fit.weights = alpha;
        fit.pvals = stats.p;
        fit.se = stats.se;
     
        %A =[ones(size(Xfull,1),1) Xfull];
        A = [Xfull]; %%% constant off
        y_hat = A*round(alpha,10);
        error = (y_hat - y).^2;
        rss = sum(error);

        % Fit partial model and calculate RSS0 for all minus subset
        for ii = 1:nelectrodes
          indices = ~(e_names == ii); % returns logical vector
          X0 = X(:,indices); % only look at subset of history not including electrode ii
          X0full = X0 * Z0;  % regressors for y_hat = X0*Z0*a0
          [a0,~,stats] = glmfit(X0full,y,'normal','constant','off');
          fit0.weights = a0;
          fit0.pvals = stats.p;
          fit0.se = stats.se;

         % A =[ones(size(X0full,1),1) X0full];
           A = [X0full];
          y_hat = A*round(a0,10) ;

          error = (y_hat - y).^2;
          rss0 = sum(error);
          
          % Compute F statistic
  
           F(electrode,ii) = ((rss0 - rss)/num_c_pts)/(rss/(nobservations-nelectrodes*num_c_pts));
          
        end
        
    end
 
    
    % Hypothesis test
    
   %  F(1:size(F,1)+1:end) = NaN;
    adj_mat = fpdf(F,num_c_pts,nobservations-nelectrodes*num_c_pts);
    
    %q = 0.1; % max number acceptable proportion of false discoveries 
  %  m = nelectrodes^2 - nelectrodes; % number of total tests performed
     m = nelectrodes^2;
    ivals = 1:m;
    sp = ivals*q/m;
    [pvals, index] = sort(adj_mat(:));
  %   pvals(isnan(pvals)) = [];
    R = find(sp>pvals'); % indices to reject null
    adj_mat = zeros(nelectrodes);
    adj_mat(index(R)) = 1; % reject H0 -> correlation
    
    
   
end

