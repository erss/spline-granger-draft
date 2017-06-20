
%%%%%%%% Three node network simulations -----------------------------------
clear all;
 close all
%%% Define model inputs ---------------------------------------------------

global s;               % tension parameter
s = 0.5;

global nsurrogates;     % number of surrogates
nsurrogates = 1000;

nelectrodes = 9; % number of electrodes
nlags = 100;  % true model order

T = 2;      % total length of recording (seconds)


f0 = 500;  % sampling frequency (Hz)
dt = 1/f0; % seconds
df = 1/T;   % frequency resolution
fNQ = f0/2; % Nyquist frequency

N = T*f0 + nlags + 100;
taxis = dt:dt:T; % time axis
noise = 0.5;
data = zeros(nelectrodes,N);
model_order = 50; % order used in model estimation

%coefficients_9N; % coefficients
load('b_50lag_sparse.mat');
b=bhat;
clear bhat;
%%% Plot true coefficients
figure;
[bb, ii] = max(abs(b),[],3);
%bb(bb<0.035)=0;
%bb(bb>=0.035)=1;
plotNetwork(bb)

%%% Simulate data ------------------------------------------
for k = nlags:length(data)-1
    data(:,k+1) = myPrediction(data(:,1:k),b);
    data(:,k+1) = data(:,k+1) + noise.*randn(size(data,1),1);
end
data = data(:,nlags+101:end);

%%% Fit standard AR to data ----------------------------------------------
%  X=data;
%  mvar_aic
tic
[ adj_standard] = build_ar( data, model_order);
standardtime  = toc;
%%% Fit spline to data ---------------------------------------------------

cntrl_pts = make_knots(model_order,floor(model_order/3));
tic
[ adj_mat] = build_ar_splines( data, model_order, cntrl_pts );
splinetime  = toc;
[ bhat, yestimate ] = estimate_coefficient_fits( data, adj_mat, model_order, cntrl_pts);
[b_est_stand, y_est_stand] = estimate_standard( data, adj_standard, model_order );

% 
for i = 1:nelectrodes
   figure;
   subplot 121
   mySpec(data(i,:),f0,'yesplot','tapers');
   subplot 122
   mySpec(yestimate(i,:),f0,'yesplot','tapers');
end
%%% Plot results ----------------------------------------------------------
figure;
subplot(3,3,[1 3])
plotchannels(taxis,data')
subplot(3,3,4)
adj_true =bb;
adj_true=b;
adj_true(adj_true~=0)=1;
adj_true=sum(adj_true,3);
adj_true(adj_true~=0)=1;
plotNetwork(adj_true)
title('True Network')

h1=subplot(3,3,5);
plotNetwork(adj_standard)
title(strcat({'Standard, '},num2str(standardtime),{' s'}))
dj = jdist(adj_true,adj_standard);

h2 = subplot(3,3,8);
xl = xlim(h2); 
xPos = xl(1) + diff(xl) / 2; 
yl = ylim(h2);
yPos = yl(2)-0.2;
t = text(xPos, yPos, sprintf('%s\n%s\n%s', 'Jaccard Dist:', num2str(dj)), 'Parent', h2);
set(t, 'HorizontalAlignment', 'center','FontSize',13);
set ( h2, 'visible', 'off')

h1 = subplot(3,3,6);
plotNetwork(adj_mat)
title('Spline Network')
title(strcat({'Spline, '},num2str(splinetime),{' s'}))
dj = jdist(adj_true,adj_mat);
h2 = subplot(3,3,9);
xl = xlim(h2); 
xPos = xl(1) + diff(xl) / 2; 
yl = ylim(h2)-0.2;
yPos = yl(2);
t = text(xPos, yPos, sprintf('%s\n%s\n%s', 'Jaccard Dist:', num2str(dj)), 'Parent', h2);
set(t, 'HorizontalAlignment', 'center','FontSize',13);
set ( h2, 'visible', 'off')


% goodness_of_fit_residuals;
% goodness_of_fit_bootstrap;
% goodness_of_fit_spectrum;

