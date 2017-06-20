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

if ~exist('data_type','var')
  data_type = 'simulation';
end

if ~exist('b_est_stand','var')
  b_est_stand = NaN(size(b));
end


X=[];
for electrode = 1:nelectrodes % run GOF on all electrodes
    if strcmp(data_type,'real') || strcmp(noise_type,'pink')
       
        data_true=data;
        [faxis, h] = mySpec(data(electrode,:),f0,'noplot','notapers' );
    else
    % Simulate data using true coefficients 
        h_sum = 0;
        for i = 1:nrealizations
            data_true = data;
             if strcmp(noise_type,'white')

                data_true = zeros(nelectrodes,N);
                for k = nlags:length(data_true)-1
                    data_true(:,k+1) = myPrediction(data_true(:,1:k),b);
                    data_true(:,k+1) = data_true(:,k+1) + noise.*randn(nelectrodes,1);
                end
                data_true= data_true(:,nlags+1:end);

% 
%             elseif strcmp(noise_type,'pink')
%                  alpha = 0.33;
%                  data_true  = make_pink_noise(alpha,nobs,dt);

             end

             y = data_true(electrode,:);   
            [faxis, h] = mySpec( y, f0,'noplot','notapers' );
            h_sum = h + h_sum;
        end
           h = h_sum/nrealizations;
    end
    % Simulate data using estimated SPLINE % STAND coefficients   
    
 
              h_sum = 0;
              h_sum_stand = 0;
        for i = 1:nrealizations
          
            
                data_hat = zeros(nelectrodes,N+model_order-nlags);
                data_stand = zeros(nelectrodes,N+model_order-nlags);
                for k = model_order:length(data_hat)-1
                    data_hat(:,k+1) = myPrediction(data_hat(:,1:k),bhat);
                    data_hat(:,k+1) = data_hat(:,k+1) + noise.*randn(nelectrodes,1);
                    
                    data_stand(:,k+1) = myPrediction(data_stand(:,1:k),b_est_stand);
                    data_stand(:,k+1) = data_stand(:,k+1) + noise.*randn(nelectrodes,1);
                    
                    
                end
               data_hat= data_hat(:,model_order+1:end);
               data_stand= data_stand(:,model_order+1:end);
                yhat = data_hat(electrode,:);   
                [faxis_hat, h_hat] = mySpec( yhat, f0,'noplot' ,'notapers'); % compute spectra
                h_sum = h_hat + h_sum;
                
                
                  ystand = data_stand(electrode,:);   
                [faxis_stand, h_stand] = mySpec( ystand, f0,'noplot' ,'notapers'); % compute spectra
                h_sum_stand = h_stand + h_sum_stand;
        end
        h_hat = h_sum/nrealizations;
        h_stand = h_sum_stand/nrealizations;
     

    % Plot spectra 
    y = data(electrode,:);           % signal 'electrode' in 'true network'
    yhat = data_hat(electrode,:);    % signal 'electrode' in 'estimated network'
    ystand = data_stand(electrode,:);
    figure;
    subplot 231
    %[faxis, h] = mySpec( y, f0 );
    plot(faxis,h);     
    xlim([0 f0/4]);
    xlabel('Frequency (Hz)','FontSize',15);
    ylabel('Averaged Power','FontSize',15);
    title('true signal','FontSize',15);

    subplot 232
    %[faxis_hat, h_hat] = mySpec( yhat, f0 );
    %title('estimated signal');
    plot(faxis_hat,h_hat);     
    xlim([0 f0/4]);
    xlabel('Frequency (Hz)','FontSize',15);
    ylabel('Averaged Power','FontSize',15);
    title('spline signal','FontSize',15);
    str1 = strcat({'Averaged spectra of electrode '},num2str(electrode),{' for '}, num2str(nrealizations),{' realizations.'});
    suptitle(str1);
    
     subplot 233
    %[faxis_hat, h_hat] = mySpec( yhat, f0 );
    %title('estimated signal');
    plot(faxis_stand,h_stand);     
    xlim([0 f0/4]);
    xlabel('Frequency (Hz)','FontSize',15);
    ylabel('Averaged Power','FontSize',15);
    title('standard signal','FontSize',15);
    str1 = strcat({'Averaged spectra of electrode '},num2str(electrode),{' for '}, num2str(nrealizations),{' realizations.'});
    suptitle(str1);

    % Compute empircal cdf
    [H, X] = ecdf(h);           
    [H1, X1] = ecdf(h_hat);
    [H2, X2] = ecdf(h_stand);


    % 
    % H=2*H;
    % H1 = 2*H1;



    subplot(2,3,[4 5 6])
    plot(X1,H1,'r','LineWidth',1.5);
    hold on
    plot(X,H,'k','LineWidth',1.5);
        plot(X2,H2,'g','LineWidth',1.5);
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
    
    
    
    R = xcov(ystand,flag); % autocovariance of estimated signal

    G = sum(R(3:end-2).^2);
    G = G/(4*pi);

    conf2 = ap*sqrt(8*pi*G/total_observations);
    UB = H2 + conf2;
    LB = H2 - conf2;
    plot(X2,UB, '--g');
    plot(X2,LB, '--g');
    axis tight
    
    
    legend('ESpline,Network','True Network','Standard Network')
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
    %str1 = strcat({'CDFs of Averaged Spectrum, '},num2str(percent_in_bounds),{' % in bds '});
    title('CDFs of Averaged Spectrum','FontSize',15);
    
    % pg 476
    gr_statistic = max(sqrt(total_observations)*abs(Qp(2,:)-Qp(4,:)))
    if (gr_statistic < ap)
        fprintf('good fit %d\n')        

    else
        fprintf('bad fit %d\n')
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
