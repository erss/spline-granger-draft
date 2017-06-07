clear all
close all
%%%% ---- Load data --------------------------------------------------


load('data_unfiltered.mat')
data = data - repmat(mean(data,2), [1,129290]); % remove mean
data_type = 'real';
N = 1000; % number of samples, 2 seconds

sz = 9e4:9e4+N-1;      % seizure
pre_sz = 2e4:2e4+N-1;  % preseizure

%signal = [2 4 6 18 22 34 36 38 57 64];

signal = [38];
data_sz= data(signal,sz);
f0 = 499.7;     % sampling frequency (Hz)

X = data(signal,sz);  %%% AIC test
amo= NaN;
mvar_aic;

model_order = 10;%morder;     %%% Choose model order
cntrl_pts = make_knots(model_order,floor(model_order/3));

%%%% ---- Define inputs --------------------------------------------------


T = N/f0;         % total length of recording (seconds)
dt = 1/f0;    % seconds

df = 1/T;      % frequency resolution
fNQ = f0/2;    % Nyquist frequency
noise= 0.25;
taxis = dt:dt:T; % time axis

global s;
s = 0.5;
global nsurrogates;
nsurrogates = 10000;

nlags = model_order;

%%% ---- Build network --------------------------------------------------

tic
[ adj_standard ] = build_ar( data_sz, model_order); % Build network using splines
standard_time  = toc;
fprintf('\nMVGC time = %d\n',standard_time);
fprintf('using model order = %d\n',model_order);
tic
[ adj_spline ] = build_ar_splines( data_sz, model_order,cntrl_pts); % Build network using splines
spline_time  = toc;
fprintf('\nspline MVGC time = %d\n',spline_time);
fprintf('using model order = %d\n',model_order);



[ bhat, yestimate ] = estimate_coefficient_fits( data_sz, adj_spline, cntrl_pts(end), cntrl_pts);
nelectrodes = length(signal);
adj_mat = adj_spline;

data = data_sz;
adj_true = adj_spline;
b= zeros(nelectrodes,nelectrodes,model_order);

figure;
subplot (2,2,[1 2])
plot(taxis,data);
subplot (2,2,3)
mySpec(data,f0,'yesplot','tapers')
title('true spectrum')
subplot (2,2,4)
mySpec(yestimate,f0,'yesplot','tapers')
title('estimated spectrum')

% 
 goodness_of_fit_residuals;
 goodness_of_fit_bootstrap;
% 




%%%% Make table of all inputs ------------------------
Sampling_Frequency = f0;
Noise_Variance = noise.^2;
T_seconds = T;
True_Model_Order = nlags;
Estimated_Order = model_order;
DW_statistic_e1 = dw(1);
% DW_statistic_e2 = dw(2);
% DW_statistic_e3 = dw(3);
pval_e1 = pval(1);
% pval_e2 = pval(2);
% pval_e3 = pval(3);
network = num2str(signal);
Tension_Parameter = s;
Tp = table(Sampling_Frequency,Noise_Variance,T_seconds,signal, True_Model_Order,...
    Estimated_Order,Tension_Parameter,moAIC,moBIC, DW_statistic_e1, pval_e1,...
    spline_time, standard_time);
% Tp = table(Sampling_Frequency,Noise_Variance,T_seconds, True_Model_Order,...
%     Estimated_Order,Tension_Parameter,moAIC,moBIC, DW_statistic_e1, pval_e1,...
%     DW_statistic_e2, pval_e2, DW_statistic_e3, pval_e3,spline_time, standard_time);
figure;
uitable('Data',Tp{:,:}','RowName',Tp.Properties.VariableNames,...
    'Units', 'Normalized', 'Position',[0,0,1,1]);


% % Save all simulation and table plots ---------------------------


h = get(0,'children');
for i=1:length(h)
    saveas(h(i), ['1N_realdata_modelorder_' num2str(model_order) '_' num2str(i)], 'jpg');
end