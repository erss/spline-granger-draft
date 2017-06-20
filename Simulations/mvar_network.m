clear all;
close all;
ntrials   = 1;     % number of trials
nobs      = 1000;   % number of observations per trial

regmode   = 'OLS';  % VAR model estimation regression mode ('OLS', 'LWR' or empty for default)
icregmode = 'LWR';  % information criteria regression mode ('OLS', 'LWR' or empty for default)

morder    = 'AIC';  % model order to use ('actual', 'AIC', 'BIC' or supplied numerical value)
momax     = 20;     % maximum model order for model order estimation

acmaxlags = 1000;   % maximum autocovariance lags (empty for automatic calculation)

tstat     = '';     % statistical test for MVGC:  'F' for Granger's F-test (default) or 'chi2' for Geweke's chi2 test
alpha     = 0.05;   % significance level for significance test
mhtc      = 'FDR';  % multiple hypothesis test correction (see routine 'significance')

fs        = 200;    % sample rate (Hz)
fres      = [];     % frequency resolution (empty for automatic calculation)

seed      = 0;      % random seed (0 for unseeded)

% Generate VAR test data (<mvgc_schema.html#3 |A3|>)
%
% _*Note:*_ This is where you would read in your own time series data; it should
% be assigned to the variable |X| (see below and <mvgchelp.html#4 Common
% variable names and data structures>).

% Seed random number generator.

rng_seed(seed);


% Get VAR coefficients for 5-node test network.

AT = var9_test; %var5_test;
nvars = size(AT,1); % number of variables
a = AT;
amo=3;
% Residuals covariance matrix.

SIGT = eye(nvars);

% Generate multi-trial VAR time series data with normally distributed residuals
% for specified coefficients and covariance matrix.

ptic('\n*** var_to_tsdata... ');
X = var_to_tsdata(AT,SIGT,nobs,ntrials);
ptoc;

%
data = X;
model_order = 40;
global s;
s = 0.5;
global nsurrogates;
nsurrogates = 100;
f0=fs;

noise=1;
N=size(data,2);
T= N/f0;
dt=1/f0;
taxis = dt:dt:T; % time axis
nlags = 3;
nelectrodes = size(data,1);
%%% Fit spline to data ---------------------------------------------------

cntrl_pts = make_knots(model_order,10);
tic
[ adj_mat] = build_ar_splines( data, model_order, cntrl_pts );
splinetime  = toc;
[ bhat, yestimate ] = estimate_coefficient_fits( data, adj_mat, model_order, cntrl_pts);
%mvar_aic;
%%
tic
[ adj_standard] = build_ar( data, model_order);
standardtime = toc;
[ b_est_stand, ~ ] = estimate_standard( data, adj_mat, model_order);

%%% Plot results ----------------------------------------------------------
% subplot(2,3,[1 3])
% for i = 1:9
%    plot(taxis,data(i,:));
%  hold on; 
% end

ylabel('Signal')
xlabel('Time (seconds)')

title('Simulated Signal','FontSize',15);
b=AT;
b(b~=0)=1;
b =sum(b,3);
b(b~=0)=1;
adj_true = sum(AT,3);
adj_true(adj_true~=0) = 1;


subplot(1,3,1)
plotNetwork(b)
axis square
title('True Network','FontSize',15)
b=AT;
subplot(1,3,2)
plotNetwork(adj_standard)
axis square

title(strcat({'Standard MVGC, '},num2str(standardtime),{' s'}),'FontSize',15)

subplot(1,3,3)
plotNetwork(adj_mat)
axis square
title(strcat({'Spline MVGC, '},num2str(splinetime),{' s'}),'FontSize',15)

% Sampling_Frequency = f0;
% Noise_Variance = noise.^2;
% T_seconds = T;
% Model_Order = nlags;
% Estimated_Order = model_order;
% Tension_Parameter = s;
% Tp = table(Sampling_Frequency,Noise_Variance,T_seconds,Model_Order,...
%     Estimated_Order,Tension_Parameter);
% figure;
% uitable('Data',Tp{:,:},'ColumnName',Tp.Properties.VariableNames,...
%     'RowName',Tp.Properties.RowNames,'Units', 'Normalized', 'Position',[0, 0, 1, 1]);
% 
% b=AT;
%%% Plot all results --------------------------------------------

% Save all simulation and table plots ---------------------------
% 
% h = get(0,'children');
% j=1;
% for i=length(h):-1:1
%     saveas(h(j), ['9N_'  num2str(i) '_summaryplot' num2str(i)], 'jpg');
%     j=j+1;
% end
% close all
% 
% % Spectral GoF --------------------------------------------------
% goodness_of_fit_spectrum;
% h = get(0,'children');
% j=1;
% for i=length(h):-1:1
%     saveas(h(j), ['9N_MVAR_e' num2str(i) '_spectrum'], 'jpg');
%     j=j+1;
% end
% close all
% 
% % Boostrap GoF --------------------------------------------------
% goodness_of_fit_bootstrap;
% h = get(0,'children');
% j=1;
% for i=length(h):-1:1
%     saveas(h(j), ['9N_MVAR_e' num2str(i) '_bootstrap'], 'jpg');
%     j=j+1;
% end
% close all
% % Residuals GoF --------------------------------------------------
% goodness_of_fit_residuals;
% h = get(0,'children');
% j=1;
% for i=length(h):-1:1
%     saveas(h(j), ['9N_MVAR_e' num2str(i) '_residuals'], 'jpg');
%     j=j+1;
% end
% close all


% goodness_of_fit_residuals;
% goodness_of_fit_bootstrap;
% goodness_of_fit_spectrum;

