
%%%%%%%% Three node network simulations -----------------------------------
clear all;
close all
%%% Define model inputs ---------------------------------------------------

global s;               % tension parameter
s = 0.5;

global nsurrogates;     % number of surrogates
nsurrogates = 10000;

nelectrodes = 9; % number of electrodes
nlags = 20;  % true model order
model_order = 40; % order used in model estimation
T = 2;      % total length of recording (seconds)


f0 = 500;  % sampling frequency (Hz)
dt = 1/f0; % seconds
df = 1/T;   % frequency resolution
fNQ = f0/2; % Nyquist frequency

N = T*f0 + nlags;
taxis = dt:dt:T; % time axis
noise = .05;
data = zeros(nelectrodes,N);

coefficients_9N;
%b(7,7,:) = 0.3*b(7,7,:);
b=0.6*b;
 b(8,:,:) = -0.7*b(8,:,:);
b(7,:,:) = .9*b(7,:,:);
 b(5,:,:) = 0.7*b(5,:,:);
 b(1,:,:) = 0.5*b(1,:,:);

%%% Simulate data ------------------------------------------
for k = nlags:length(data)-1;
    data(:,k+1) = myPrediction(data(:,1:k),b);
    data(:,k+1) = data(:,k+1) + noise.*randn(size(data,1),1);
end
data = data(:,nlags+1:end);

subplot(2,3,[1 3])
for i = 1:nelectrodes
plot(dt:dt:T,data(i,:));
hold on;
end
figure; plotchannels(data')

ylabel('Signal')
xlabel('Time (seconds)')
title('Simulated Signal','FontSize',15);
%%

% subplot(2,2,3)
% mySpec(data(1,:),f0);

%%% Fit standard AR to data ----------------------------------------------
X=data;
mvar_aic
tic
[ adj_standard] = build_ar( data, morder);
standardtime  = toc;
%%% Fit spline to data ---------------------------------------------------

cntrl_pts = make_knots(model_order,floor(model_order/3));
tic
[ adj_mat] = build_ar_splines( data, model_order, cntrl_pts );
splinetime  = toc;
[ bhat, yestimate ] = estimate_coefficient_fits( data, adj_mat, model_order, cntrl_pts);

for i = 1:nelectrodes
   figure;
   subplot 121
   mySpec(data(i,:),f0,'yesplot','tapers');
   subplot 122
   mySpec(yestimate(i,:),f0,'yesplot','tapers');
end

%%% Plot results ----------------------------------------------------------

subplot(2,3,[1 3])
plotchannels(data')
subplot(2,3,4)
adj_true=b;
adj_true(adj_true~=0)=1;
adj_true=sum(adj_true,3);
adj_true(adj_true~=0)=1;
plotNetwork(adj_true)
title('True Network')

subplot(2,3,5)
plotNetwork(adj_standard)
title(strcat({'Standard, '},num2str(standardtime),{' s'}))


subplot(2,3,6)
plotNetwork(adj_mat)
title('Spline Network')
title(strcat({'Spline, '},num2str(splinetime),{' s'}))



% goodness_of_fit_residuals;
% goodness_of_fit_bootstrap;
% goodness_of_fit_spectrum;


X=data;
mvar_aic;