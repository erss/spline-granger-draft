%%% Compare spectrum of true signal and of estimated signal
%%%%%% NOTE highly dependent on noise used to generate model
% Define model inputs
electrode = 1; % which electrode to run GOF on
nrealizations = 1; % number of realizations for each process


% Simulate data using true coefficients 
h_sum = 0;
for i = 1:nrealizations
    data_true = zeros(nelectrodes,N);
    for k = nlags:length(data_true)-1;
        data_true(:,k+1) = myPrediction(data_true(:,1:k),b);
        data_true(:,k+1) = data_true(:,k+1) + noise.*randn(nelectrodes,1);
    end
    
    data_true= data_true(:,41:end);
    
     y = data_true(electrode,:);   
    [faxis, h] = mySpec( y, f0,0 );
    h_sum = h + h_sum;
end
   h = h_sum/nrealizations;
   
% Simulate data using estimated coefficients   
   h_sum = 0;
for i = 1:nrealizations
        data_hat = zeros(nelectrodes,N);
    for k = model_order:length(data_hat)-1;
        data_hat(:,k+1) = myPrediction(data_hat(:,1:k),bhat);
        data_hat(:,k+1) = data_hat(:,k+1) + noise.*randn(nelectrodes,1);
    end
    data_hat= data_hat(:,nlags+1:end);
    yhat = data_hat(electrode,:);   
    [faxis_hat, h_hat] = mySpec( yhat, f0,0 ); % compute spectra
    h_sum = h_hat + h_sum;
end
h_hat = h_sum/nrealizations;


%%% Simulate independent network ---------------------------------------

a1 = 0.07*[hann(20)', -0.5*ones(20,1)']';   %AR coeffictients for signal 1
a2 = 0.03*[-0.5*ones(20,1)', hann(20)']';   %                  ...signal 2
a3 = -.3*ones(size(a1));                    %                  ...signal 3

z = zeros(3,3*nlags);
z(1,1:40) = a1;
z(1,41:80) = a2;
z(2,1:40) = a1;
z(2,41:80) = a2;
z(3,81:120) = a3;


h_sum = 0;
for i = 1:nrealizations
    data_z = zeros(3,N);
    
    for k = nlags:length(data_z)-1;
        data_z(:,k+1) = myPrediction(data_z(:,1:k),z);
        data_z(:,k+1) = data_z(:,k+1) + noise.*randn(3,1);
    end
    
    data_z= data_z(:,41:end);
    zhat = data_z(electrode,:);    
    [faxis_z, h_z] = mySpec( zhat, f0,0 );
    h_sum = h_z + h_sum;
end
h_z = h_sum/nrealizations;


% Plot spectra 
y = data(electrode,:);           % signal 1 in 'true network'
yhat = data_hat(electrode,:);    % signal 1 in 'estimated network'
zhat = data_z(electrode,:);      % signal 1 in 'independent network'

figure;
subplot 121
%[faxis, h] = mySpec( y, f0 );
plot(faxis,h);     
xlim([0 f0/4]);
xlabel('Frequency (Hz)','FontSize',15);
ylabel('Power','FontSize',15);
title('true signal','FontSize',15);

subplot 122
%[faxis_hat, h_hat] = mySpec( yhat, f0 );
%title('estimated signal');
plot(faxis_hat,h_hat);     
xlim([0 f0/4]);
xlabel('Frequency (Hz)','FontSize',15);
ylabel('Power','FontSize',15);
title('true signal','FontSize',15);



% Compute empircal cdf
[H, X] = ecdf(h);           
[H1, X1] = ecdf(h_hat);


% 
% H=2*H;
% H1 = 2*H1;



figure();
plot(X1,H1,'r','LineWidth',1.5);
hold on
plot(X,H,'k','LineWidth',1.5);


%legend('Estimated Signal','True Signal')
title('CDFs of Spectrum');

% Compute confidence bounds for estimated signal (Priestley p 478)

a = 2.2414; % for 95% confidence bounds
N = size(data_true,2); %check! number of observations from which H is computed ? length of signal ?

flag = 'biased'; % divide by 1/N

R = xcov(yhat,flag); % autocovariance of estimated signal

G = sum(R(3:end-2).^2);
G = G/(4*pi);

conf1 = a*sqrt(8*pi*G/N);

plot(X1,H1 + conf1, '--r');
plot(X1,H1 - conf1, '--r');


