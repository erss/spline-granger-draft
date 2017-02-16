%%% Spectral confidence bounds
close all; clear all;
T = 5;         % total length of recording (seconds)
dt = 0.001;    % seconds

f0 = 1/dt;     % sampling frequency (Hz)
N1 = T*f0;     % number of samples needed
df = 1/T;      % frequency resolution
fNQ = f0/2;    % Nyquist frequency
noise=.25;
taxis = dt:dt:T; % time axis

electrode = 1; % which electrode to run GOF on
nrealizations = 1; % number of realizations for each process

%%% simulate 'true' network --------------------------------------------------
a1 = 0.07*[hann(20)', -0.5*ones(20,1)']';   %AR coefficients for signal 1
a2 = 0.05*[-0.5*ones(20,1)', hann(20)']';   %                  ...signal 2
a3 = -.3*ones(size(a1));                    %                  ...signal 3

L = length(a1);                             % Number of AR terms.
N = N1+L;                                   % Number of time steps.
           
nlags = 40;                               % Define order of AR model
                                           % needs to be larger than true
                                           % order                                      
b = zeros(3,3*nlags);
b(1,1:40) = a2;
b(1,41:80) = a3;
b(2,41:80) = a3;
b(3,81:120) = a1;

h_sum = 0;
for i = 1:nrealizations
    data = zeros(3,N);
    for k = nlags:length(data)-1;
        data(:,k+1) = myPrediction(data(:,1:k),b,nlags);
        data(:,k+1) = data(:,k+1) + noise.*randn(3,1);
    end
    
    data= data(:,41:end);
   y = data(electrode,:);   
    [faxis, h] = mySpec( y, f0,0 );
    h_sum = h + h_sum;
end
   h = h_sum/nrealizations;
figure;
subplot 311
 plot(taxis,data(1,:));
 hold on
 plot(taxis,data(2,:));
 plot(taxis,data(3,:));
ylabel('Signal')
xlabel('Time (seconds)')
legend('x1','x2','x3')
title('true network','FontSize',15);

%%% simulate 'estimated' network------------------------------------
% Estimate network using splines
[ adj_mat ] = build_ar_splines( data, nlags); % Build network using splines

%%Get coefficient estimates and signal estimates

flag = 1; % use splines to estimate
[bhat, yhat] = estimate_coef(data,adj_mat,nlags,flag);

h_sum = 0;
for i = 1:nrealizations
        data_hat = zeros(3,N);
    for k = nlags:length(data_hat)-1;
        data_hat(:,k+1) = myPrediction(data_hat(:,1:k),bhat,nlags);
        data_hat(:,k+1) = data_hat(:,k+1) + noise.*randn(3,1);
    end
    data_hat= data_hat(:,41:end);
    yhat = data_hat(electrode,:);   
    [faxis, h_hat] = mySpec( yhat, f0,0 ); % compute spectra
    h_sum = h_hat + h_sum;
end
h_hat = h_sum/nrealizations;

subplot 312
 plot(taxis,data_hat(1,:));
 hold on
 plot(taxis,data_hat(2,:));
 plot(taxis,data_hat(3,:));
ylabel('Signal')
xlabel('Time (seconds)')
legend('x1','x2','x3')
title('estimated network','FontSize',15);


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
        data_z(:,k+1) = myPrediction(data_z(:,1:k),z,nlags);
        data_z(:,k+1) = data_z(:,k+1) + noise.*randn(3,1);
    end
    
    data_z= data_z(:,41:end);
    zhat = data_z(electrode,:);    
    [faxis, h_z] = mySpec( zhat, f0,0 );
    h_sum = h_z + h_sum;
end
h_z = h_sum/nrealizations;





subplot 313
 plot(taxis,data_z(1,:));
 hold on
 plot(taxis,data_z(2,:));
 plot(taxis,data_z(3,:));
ylabel('Signal')
xlabel('Time (seconds)')
legend('x1','x2','x3')
title('independent network','FontSize',15);



%%% Construct goodness-of-fit -------------------------------------------

% Plot spectra 
y = data(electrode,:);           % signal 1 in 'true network'
yhat = data_hat(electrode,:);    % signal 1 in 'estimated network'
zhat = data_z(electrode,:);      % signal 1 in 'independent network'

figure;
subplot 131
%[faxis, h] = mySpec( y, f0 );
plot(faxis,h);     
xlim([0 f0/4]);
xlabel('Frequency (Hz)','FontSize',15);
ylabel('Power','FontSize',15);
title('true signal','FontSize',15);

subplot 132
%[faxis_hat, h_hat] = mySpec( yhat, f0 );
%title('estimated signal');
plot(faxis,h_hat);     
xlim([0 f0/4]);
xlabel('Frequency (Hz)','FontSize',15);
ylabel('Power','FontSize',15);
title('true signal','FontSize',15);

subplot 133
%[faxis_z, h_z] = mySpec( zhat, f0 ,1);
%title('independent signal');
plot(faxis,h_z);     
xlim([0 f0/4]);
xlabel('Frequency (Hz)','FontSize',15);
ylabel('Power','FontSize',15);
title('true signal','FontSize',15);

% compute and plot cumulative distributions of spectra
[H, X] = ecdf(h);           
[H1, X1] = ecdf(h_hat);
[H2, X2] = ecdf(h_z);
% 
% H=2*H;
% H1 = 2*H1;
% H2 = 2*H2;


figure();
plot(X1,H1,'r','LineWidth',1.5);
hold on
plot(X2,H2,'g','LineWidth',1.5);
plot(X,H,'k','LineWidth',1.5);
legend('Estimated Signal','Independent Signal','True Signal')
title('CDFs of Spectrum');

% Compute confidence bounds for estimated signal (Priestley p 478)

a = 2.2414; % for 95% confidence bounds
N = N1; % number of observations from which H is computed ? length of signal ?

flag = 'biased'; % divide by 1/N

R = xcov(yhat,flag); % autocovariance of estimated signal

G = sum(R(3:end-2).^2);
G = G/(4*pi);

conf1 = a*sqrt(8*pi*G/N);

plot(X1,H1 + conf1, '--r');
plot(X1,H1 - conf1, '--r');


% Compute confidence bounds for independent signal (Priestley p 478)

R = xcov(zhat,flag); % autocovariance of independent signal
G = sum(R(3:end-2).^2); % 3:9997, M = 3
G = G/(4*pi);

conf2 = a*sqrt(8*pi*G/N);

plot(X2,H2 + conf2, '--g');
plot(X2,H2 - conf2, '--g');

