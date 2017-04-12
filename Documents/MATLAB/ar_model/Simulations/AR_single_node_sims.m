%%%%%%%% Single node network simulations ----------------------------------
clear all;
%%% Define model inputs ---------------------------------------------------
nobs = 1000;   % number of observations per trial
sgnl = 1;
nelectrodes = 1;

f0 = 200;  % sampling frequency (Hz)
T = nobs/f0;      % total length of recording (seconds)
dt = 1/f0; % seconds

N = T*f0 +2;   % number of samples needed +2
df = 1/T;   % frequency resolution
fNQ = f0/2; % Nyquist frequency

taxis = dt:dt:T; % time axis
noise = 0.7;
data = zeros(1,N);
model_order = 20; % order used in model estimation

%%%% WHITE NOISE ---------------------------------------------------------
%%% Sim 1 - HIGH FREQUENCY -----------------------------------------------
b = zeros(1,1,2);
b(1,1,:)= [0.9 -0.8];
nlags = size(b,3);
%%% Note AIC/BIC indicate 2 is the best order

% %%% Sim 2 - LOW FREQUENCY I -----------------------------------------------
% b = zeros(1,1,1);
% b(1,1,:)= [0.9];
% nlags = size(b,3);
% 
% %%% Sim 3 - LOW FREQUENCY II ----------------------------------------------
% b = zeros(1,1,2);
% b(1,1,:)= [0.3 0.3];
% nlags = size(b,3);
%  
% %%% Sim 3 - LOW FREQUENCY III ---------------------------------------------
% 
% b = zeros(1,1,2);
% b(1,1,:)= [0.9 -0.1];
% nlags = size(b,3);


%%% Generate white noise data from above coefficients
for k = nlags:length(data)-1;
   data(:,k+1) = myPrediction(data(:,1:k),b);
   data(:,k+1) = data(:,k+1) + noise.*randn(size(data,1),1);
      
end

%%%% PINK NOISE ----------------------------------------------------------
%  alpha = 0.33;
%  data  = make_pink_noise(alpha,nobs,dt);

subplot(3,2,[1 2])
 plot(data(1,:));

ylabel('Signal')
xlabel('Time (seconds)')
legend('x1')
title('Simulated Signal','FontSize',15);


subplot(3,2,3)
mySpec(data(1,:),f0);

%%% Fit spline to data ---------------------------------------------------

cntrl_pts = make_knots(model_order,10);
[ adj_mat] = build_ar_splines( data, model_order, cntrl_pts );
[bhat, yhat] = estimate_coefficient_fits( data, adj_mat, model_order,cntrl_pts );


subplot(3,2,[5 6])
plot(squeeze(bhat(1,1,:)),'LineWidth',1.5)
hold on
plot(cntrl_pts(2:end),squeeze(bhat(1,1,cntrl_pts(2:end))),'o')
title('Estimated Coefficients','FontSize',15);
% figure;
% plot(data);
% hold on
% plot(yhat,'--r');


subplot(3,2,4)
mySpec(yhat,f0);
title('Estimated signal spectrogram','FontSize',15);


%%% Determine what AIC thinks is best order
a=b;
mvar_aic;
