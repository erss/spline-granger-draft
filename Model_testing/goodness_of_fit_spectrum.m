%%% Compare spectrum of true signal and of estimated signal
%%%%%% NOTE highly dependent on noise used to generate model
% Define model inputs

if nelectrodes > 3
    nrealizations = 10;
else
    nrealizations = 10;  % number of realizations for each process
end

if ~exist('noise_type','var')
  noise_type = 'white';
end

for electrode = 1:nelectrodes % run GOF on all electrodes
    
    % Simulate data using true coefficients 
    h_sum = 0;
    for i = 1:nrealizations
        
         if strcmp(noise_type,'white')

            data_true = zeros(nelectrodes,N);
            for k = nlags:length(data_true)-1;
                data_true(:,k+1) = myPrediction(data_true(:,1:k),b);
                data_true(:,k+1) = data_true(:,k+1) + noise.*randn(nelectrodes,1);
            end
            data_true= data_true(:,nlags+1:end);

             
        elseif strcmp(noise_type,'pink')
             alpha = 0.33;
             data_true  = make_pink_noise(alpha,nobs,dt);

         end
        
         y = data_true(electrode,:);   
        [faxis, h] = mySpec( y, f0,0 );
        h_sum = h + h_sum;
    end
       h = h_sum/nrealizations;

    % Simulate data using estimated coefficients   
    
 
              h_sum = 0;
        for i = 1:nrealizations
          
            
                data_hat = zeros(nelectrodes,N+model_order-nlags);
                
                for k = model_order:length(data_hat)-1;
                    data_hat(:,k+1) = myPrediction(data_hat(:,1:k),bhat);
                    data_hat(:,k+1) = data_hat(:,k+1) + noise.*randn(nelectrodes,1);
                end
               data_hat= data_hat(:,model_order+1:end);
      
                yhat = data_hat(electrode,:);   
                [faxis_hat, h_hat] = mySpec( yhat, f0,0 ); % compute spectra
                h_sum = h_hat + h_sum;
        end
        h_hat = h_sum/nrealizations;
  
     

    % Plot spectra 
    y = data(electrode,:);           % signal 'electrode' in 'true network'
    yhat = data_hat(electrode,:);    % signal 'electrode' in 'estimated network'

    figure;
    subplot 221
    %[faxis, h] = mySpec( y, f0 );
    plot(faxis,h);     
    xlim([0 f0/4]);
    xlabel('Frequency (Hz)','FontSize',15);
    ylabel('Averaged Power','FontSize',15);
    title('true signal','FontSize',15);

    subplot 222
    %[faxis_hat, h_hat] = mySpec( yhat, f0 );
    %title('estimated signal');
    plot(faxis_hat,h_hat);     
    xlim([0 f0/4]);
    xlabel('Frequency (Hz)','FontSize',15);
    ylabel('Averaged Power','FontSize',15);
    title('estimated signal','FontSize',15);
    str1 = strcat({'Averaged spectra of electrode '},num2str(electrode),{' for '}, num2str(nrealizations),{' realizations.'});
    suptitle(str1);

    % Compute empircal cdf
    [H, X] = ecdf(h);           
    [H1, X1] = ecdf(h_hat);


    % 
    % H=2*H;
    % H1 = 2*H1;



    subplot(2,2,[3 4])
    plot(X1,H1,'r','LineWidth',1.5);
    hold on
    plot(X,H,'k','LineWidth',1.5);
    % Compute confidence bounds for estimated signal (Priestley p 478)

    ap = 2.2414; % for 95% confidence bounds
    total_observations = size(data_true,2); %check! number of observations from which H is computed ? length of signal ?

    flag = 'biased'; % divide by 1/N

    R = xcov(yhat,flag); % autocovariance of estimated signal

    G = sum(R(3:end-2).^2);
    G = G/(4*pi);

    conf1 = ap*sqrt(8*pi*G/total_observations);
    UB = H1 + conf1;
    LB = H1 - conf1;
    plot(X1,UB, '--r');
    plot(X1,LB, '--r');
    axis tight
    legend('Estimated Network','True Network')
    %%% Compute amount of time in confidence bounds

    q1 = [X' X1'];
    Q = NaN(5,length(q1));
    Q(1,:) = q1;
    Q(2,:) = interp1q(X,H,Q(1,:)')';
    Q(3,:) = interp1q(X1,LB,Q(1,:)')';
    Q(4,:) = interp1q(X1,H1,Q(1,:)')';
    Q(5,:) = interp1q(X1,UB,Q(1,:)')'';
    
    Q = sortrows(Q',1)';
    Qp=Q';
    Qp = Qp(~any(isnan(Qp),2),:)';
    tol = 0;%0.05;
    percent_in_bounds = length(find( Qp(2,:)>=Qp(3,:)-tol & Qp(2,:) <= Qp(5,:)+tol ));
    percent_in_bounds = 100*percent_in_bounds/size(Qp,2)
    str1 = strcat({'CDFs of Averaged Spectrum, '},num2str(percent_in_bounds),{' % in bds '});
    title(str1,'FontSize',15);
    
    % pg 476
    gr_statistic = max(sqrt(total_observations)*abs(Qp(2,:)-Qp(4,:)))
    if (gr_statistic < ap)
        fprintf('good fit')
    else
        fprintf('bad fit')
    end
    dim = [0.5 0.05 0.2 0.3];
str = ['Test statistic: good fit if ' num2str(gr_statistic) '<' num2str(ap)];
annotation('textbox',dim,'String',str,'FitBoxToText','on');

   % [ii iii] = max(sqrt(total_observations)*abs(Q(2,:)-Q(4,:)))
   % plot(Q(1,iii),Q(2,iii),'o','MarkerSize',30) 
   
   
%     figure;
%     plot(X1,H1,'r','LineWidth',1.5);
%     hold on;
%     plot(X,H,'k','LineWidth',1.5);
%     plot(X1,UB, '--r');
%     plot(X1,LB, '--r');
%     
%     plot(Q(1,:),Q(2,:),'b')
%     hold on
%     plot(Q(1,:),Q(3,:),'--g')
%     plot(Q(1,:),Q(4,:),'g')
%     plot(Q(1,:),Q(5,:),'--g')
    
% figure;
% myKS(h,h_hat)
    
end
