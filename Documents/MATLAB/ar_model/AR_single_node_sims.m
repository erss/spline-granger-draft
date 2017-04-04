%%%%%%%% simulate data --------------------------------------------------
nobs = 1000;   % number of observations per trial
sgnl = 1;


% parameters
f0 = 200;  % sampling frequency (Hz)
T = nobs/f0;      % total length of recording (seconds)
dt = 1/f0; % seconds

N = T*f0;   % number of samples needed
df = 1/T;   % frequency resolution
fNQ = f0/2; % Nyquist frequency

taxis = dt:dt:T; % time axis
noise = 0.7;
data = zeros(1,N);

% % HIGH FREQUENCY SCENARIO 

b= [0.9 -0.8 .1 .1 .1];
nlags = length(b);
% % % 
% % LOW FREQUENCY SCENARIO 
% b= [0.9];
% nlags = length(b);
% 
% %LOW FREQUENCY SCENARIO 
% b= [0.3 0.3];
% nlags = length(b);
% 
% % LOW FREQUENCY SCENARIO 

% b= [0.9 -0.1];
% nlags = length(b);
 
for k = nlags:length(data)-1;
    data(:,k+1) = myPrediction(data(:,1:k),b,nlags);
    data(:,k+1) = data(:,k+1) + noise.*randn;
end
subplot(3,2,[1 2])
 plot(data(1,:));

ylabel('Signal')
xlabel('Time (seconds)')
legend('x1')
title('Simulated Signal','FontSize',15);


subplot(3,2,3)
mySpec(data(1,:),f0);

%%% Fit spline to data ---------------------------------------------------
nlags =100;
flag = 1; % use splines
cntrl_pts = make_knots(nlags,10);
%cntrl_pts = 0:10:nlags;
[ adj_mat] = build_ar_splines( data, nlags, cntrl_pts );
[ bhat, yhat ] = estimate_coef( data, adj_mat, nlags, flag,cntrl_pts);


subplot(3,2,[5 6])
plot(bhat,'LineWidth',1.5)
hold on
plot(cntrl_pts(2:end),bhat(cntrl_pts(2:end)),'o')
title('Estimated Coefficients','FontSize',15);
% figure;
% plot(data);
% hold on
% plot(yhat,'--r');


subplot(3,2,4)
mySpec(yhat,f0);
title('Estimated signal spectrogram','FontSize',15);
