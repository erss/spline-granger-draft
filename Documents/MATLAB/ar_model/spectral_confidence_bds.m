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

    data = zeros(3,N);
for k = nlags:length(data)-1;
    data(:,k+1) = myPrediction(data(:,1:k),b,nlags);
    data(:,k+1) = data(:,k+1) + noise.*randn(3,1);
end
data= data(:,41:end);
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


    data_hat = zeros(3,N);
for k = nlags:length(data_hat)-1;
    data_hat(:,k+1) = myPrediction(data_hat(:,1:k),bhat,nlags);
    data_hat(:,k+1) = data_hat(:,k+1) + noise.*randn(3,1);
end
data_hat= data_hat(:,41:end);
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

%%% Simulate signals --------------------------------------------------
a1 = 0.07*[hann(20)', -0.5*ones(20,1)']';   %AR coeffictients for signal 1
a2 = 0.03*[-0.5*ones(20,1)', hann(20)']';   %                  ...signal 2
a3 = -.3*ones(size(a1));                    %                  ...signal 3

z = zeros(3,3*nlags);
z(1,1:40) = a1;
z(1,41:80) = a2;
z(2,1:40) = a1;
z(2,41:80) = a2;
z(3,81:120) = a3;

    data_z = zeros(3,N);
for k = nlags:length(data_z)-1;
    data_z(:,k+1) = myPrediction(data_z(:,1:k),z,nlags);
    data_z(:,k+1) = data_z(:,k+1) + noise.*randn(3,1);
end
data_z= data_z(:,41:end);

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

% Compute spectra 
y = data(1,:);           % signal 1 in 'true network'
yhat = data_hat(1,:);    % signal 1 in 'estimated network'
zhat = data_z(1,:);      % signal 1 in 'independent network'

figure;
subplot 131
[faxis, h] = mySpec( y, f0 );
title('true signal');
subplot 132
[faxis_hat, h_hat] = mySpec( yhat, f0 );
title('estimated signal');
subplot 133
[faxis_z, h_z] = mySpec( zhat, f0 );
title('independent signal');

% compute and plot cumulative distributions of spectra
[H, X] = ecdf(h);
[H1, X1] = ecdf(h_hat);
[H2, X2] = ecdf(h_z);


figure();
plot(X,H,'k');
hold on
plot(X1,H1,'r');
plot(X2,H2,'g');
legend('True Signal','Estimated Signal','Independent Signal')
title('CDFs of Spectrum');

% Compute confidence bounds

a = 2.2414; % for 95% confidence bounds
N = N1; % number of observations from which H is computed ? length of signal ?

flag = 'biased'; % dive by 1/N

R = xcorr(y,flag); % autocorrelation of true signal

G = sum(R(2:end-2).^2);
G = G/(4*pi);

conf = a*sqrt(8*pi*G/N);

plot(X1,H1 + conf, '--r');
plot(X1,H1 - conf, '--r');

