function bhat= ar_gof( adj_spline,adj_stand, data, nlags, true_coeffs,taxis, which_electrode)
% AR_GOF  builds AR-Spline models of data for given adjacencey matrix and tests 
%         goodness of fit for one of the electrodes.
% 
% INPUTS:
%  ad_spline = the adjacency matrix for spline fit
%  adj_stand = the adjacency matrix for standard AR fit
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
% 


%%% Initialize variables
    nelectrodes = size(data,1);  % number electrodes
    y = data(which_electrode,nlags+1:end);   % true signal
    

%%% Get coefficient estimates and signal estimates for spline model-----
    flag =1; % spline model 
    [bhat, yhat] = estimate_coef(data,adj_spline,nlags,flag);
    yhat = yhat(which_electrode,:); % only compare which_electrode
    bhat = bhat(which_electrode,:); %  
    
    flag =0; % standard AR model
    [bhat0, yhat0] = estimate_coef(data,adj_stand,nlags,flag);
    yhat0 = yhat0(which_electrode,:); % only compare which_electrode
    bhat0 = bhat0(which_electrode,:); %  

    

%%% Plot coefficient fits -----------------------------------
 

nsurrogates = 10000; % number of surrogates

[~, bounds]= myBootstrap(data,adj_spline,nlags,which_electrode,nsurrogates);


%%% plot estimated coefficients and confidence bounds --------------------



   figure;
    for k = 1:nelectrodes

       subplot(1,nelectrodes,k)
       plot(taxis(1:nlags),true_coeffs{which_electrode,k},'k','LineWidth',1.5)
       hold on;
       plot (taxis(1:nlags),bhat((k-1)*nlags+1:k*nlags),'r','LineWidth',1.5)
       plot (taxis(1:nlags),bhat0((k-1)*nlags+1:k*nlags),'g','LineWidth',1.5)
       plot(taxis(1:nlags),(bounds(1,(k-1)*nlags+1:k*nlags)),'--r');
       plot(taxis(1:nlags),(bounds(2,(k-1)*nlags+1:k*nlags)),'--r');

       hold off
       xlabel('Time Lag (seconds)','FontSize', 15)
       
       if k == nelectrodes
           g = legend('True','Spline','Standard');
           set(g,'FontSize',15)
       end
    end

    
    
    

%%% Plot True & Estimated signal fits ---------------------------
    figure;
    
    subplot 311
    % true and estimated signal
    plot(taxis(nlags+1:end),y,'k','LineWidth',1)          %% plot true signal
    hold on;
    plot(taxis(nlags+1:end),yhat,'r','LineWidth',2)     %% plot estimated signal
    plot(taxis(nlags+1:end),yhat0,'g','LineWidth',1)
    h =legend('True Signal', 'Spline Model','Standard Model');
    set(h,'FontSize',12)
    title(strcat({'True and Estimated Signal for Electrode '},num2str(which_electrode)),'FontSize',15);        

    subplot 312
        plot(taxis(nlags+1:end),y,'k','LineWidth',1)          %% plot true signal
    hold on;
    plot(taxis(nlags+1:end),yhat0,'g','LineWidth',1)
    h =legend('True Signal','Standard Model');
    set(h,'FontSize',12)
    title(strcat({'True and Estimated Signal for Electrode '},num2str(which_electrode)),'FontSize',15);        

    subplot 313
        plot(taxis(nlags+1:end),y,'k','LineWidth',1)          %% plot true signal
    hold on;
    plot(taxis(nlags+1:end),yhat,'r','LineWidth',1)     %% plot estimated signal
    h =legend('True Signal', 'Spline Model');
    set(h,'FontSize',12)
    title(strcat({'True and Estimated Signal for Electrode '},num2str(which_electrode)),'FontSize',15);        

    
 %%% Plot Spectra ---------------------------
   figure; 
   f0 = 1/(taxis(2)-taxis(1));
   subplot 131
   mySpec(y,f0)
   title('True signal','FontSize',12);
   
   subplot 132
    mySpec(yhat,f0)
   title('Spline model','FontSize',12);
   
   subplot 133
    mySpec(yhat0,f0)
   title('Standard model','FontSize',12);

    
    
%     subplot 312
%     % cumulative residual process for signal estimate
%     plot(taxis(nlags+1:end),cumsum(y-yhat),'.','MarkerSize',10); hold on;    %% plot residuals
%     plot([taxis(nlags+1) taxis(end)],[0 0],'g','LineWidth',1.5)
%     title('Cumulative Residual Process for Signal Estimate','FontSize',15);        
%     xlabel('Time (seconds)','FontSize',15);
%     ylabel('Cumulative Residual','FontSize',15);
% 
%     subplot 313
%     % autocorrelation of residuals 
%     autocorr(y-yhat)

%%%  KS plot for signal fit of which_electrode----------------------------

%    figure;
%    myKS(y,yhat);

   
end

