function bhat= ar_gof( adj_mat, data, nlags, true_coeffs,taxis, which_electrode)
% AR_GOF  builds AR-Spline models of data for given adjacencey matrix and tests 
%         goodness of fit for one of the electrodes.
% 
% INPUTS:
%  ad_mat = the adjacency matrix
%  data   = A matrix of electode data with dimensions electrodes x
%           time
%  nlags  = The number of lags used as used for predictor variables
%  true_coeffs = the true coeffecients that generate the AR model
%  taxis  = time axis
%  which_electrode = the electrode whose model fit will be test
%
% OUTPUTS:
%   bhat = model coefficients
% 
% FIGURES:
%   x. KS plot for which_electrode reconstructed signal.
%   x. Plots of cumulative residual process and residual
%      autocorrelation plot for which_electrode.
%   x. Plots of estimated coefficients and true coefficients, along with
%      residuals and autocorrelation plot for which_electrode.
% 

%%% Initialize variables & outputs
    nelectrodes = size(data,1);            % number electrodes
    nobservations = length(data(1,nlags+1:end)); % number of observations


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


        %%% Plot True & Estimated signal fits ---------------------------
        if ii == which_electrode 
            %%% Plot coefficient fits -----------------------------------
        f = figure;    
        j =1;
        for k = 1:nelectrodes
           subplot(1,nelectrodes,k)
           plot(taxis(1:nlags),true_coeffs{ii,k},'LineWidth',1.5)
          % title('Coefficient Estimates','FontSize',15);
           hold on
           if preds(k)
           plot (taxis(1:nlags),bhat(((j-1)*nlags + 1): j*nlags),'--r','LineWidth',1.5)
         
           j = j+1;
           else
               plot(taxis(1:nlags),zeros(1,nlags),'--r','LineWidth',1.5)
          
           end
           hold off
           xlabel('Time Lag (seconds)','FontSize', 15)
        end
       
    
        figure;
        [yhat ylo yhi] = glmval(alpha,Xfull,'identity',stats); % Get signal estimate.
 %       yhat = yhat + noise.*randn(size(yhat));
        subplot 311
        
        plot(taxis(nlags+1:end),y,'k','LineWidth',1.5)          %% plot true signal
        hold on;
        plot(taxis(nlags+1:end),yhat,'r','LineWidth',1.5)       %% plot estimated signal
     %  plot(yhat+yhi,'--r','LineWidth',2) %% error bars
     %  plot(yhat-ylo,'--r','LineWidth',2) %%
        h =legend('True Signal', 'Estimated Signal');
        set(h,'FontSize',12)
        title(strcat({'True and Estimated Signal for Electrode '},num2str(ii)),'FontSize',15);        

    
        
        subplot 312
        plot(taxis(nlags+1:end),cumsum(y-yhat),'.','MarkerSize',10); hold on;    %% plot residuals
        plot([taxis(nlags+1) taxis(end)],[0 0],'g','LineWidth',1.5)
        title('Cumulative Residual Process for Signal Estimate','FontSize',15);        
        xlabel('Time (seconds)','FontSize',15);
        ylabel('Cumulative Residual','FontSize',15);
        
        subplot 313
        autocorr(y-yhat)
          
%%%  KS plot for signal fit of which_electrode-------------------------------------------
      
      
%          figure; suptitle('KS Plots for Signal Fits')
%         % p = get(gcf,'Number');
%       
% 
%      [emp_cdf, x1, L, U] = ecdf(yhat);
%      [theo_cdf x2, L2, U2] = ecdf(y);       
% 
%        n = length(y);
%         % KS plot
%         
%         plot(x1,emp_cdf,'k'); hold on; 
%         plot(x2,theo_cdf,'r'); 
%         plot(x1,L,'--k','LineWidth',1);
%         plot(x1,U,'--k','LineWidth',1);
%         xlabel(strcat({'Electrode '},num2str(ii)),'FontSize',12);
%       
%        axis tight
%        
%      
%        h= legend('Empirical CDF','Theoretical CDF', 'Lower CI', 'Upper CI');
%        set(h,'FontSize',8);
%        set(h,'Location','southeast')
       
       
       figure;
       
       myKS(y,yhat);
        
      end

end

end

