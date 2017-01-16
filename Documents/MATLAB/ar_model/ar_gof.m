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
% % 
% 
% %%% Initialize variables & outputs
    nelectrodes = size(data,1);            % number electrodes
    ii = which_electrode;


    % Generate observations for given y
        y = data(ii,nlags+1:end);   
    


[bhat, yhat] = estimate_coef(data,adj_mat,nlags);
 yhat = yhat(which_electrode,:);
 bhat = bhat(which_electrode,:);
        %%% Plot True & Estimated signal fits ---------------------------
      
            %%% Plot coefficient fits -----------------------------------
        figure;    
        
        for k = 1:nelectrodes

           subplot(1,nelectrodes,k)
           plot(taxis(1:nlags),true_coeffs{ii,k},'LineWidth',1.5)
           hold on;
           plot (taxis(1:nlags),bhat((k-1)*nlags+1:k*nlags),'--r','LineWidth',1.5)
        
           hold off
           xlabel('Time Lag (seconds)','FontSize', 15)
        end
          


        figure;
        subplot 311
        
        plot(taxis(nlags+1:end),y,'k','LineWidth',1.5)          %% plot true signal
        hold on;
        plot(taxis(nlags+1:end),yhat,'r','LineWidth',1.5)       %% plot estimated signal
    
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
       
       figure;
       
       myKS(y,yhat);
        
      

end

