%%%%%%%% Single node network simulations ----------------------------------
clear all;
%close all;


%%% Define model inputs ---------------------------------------------------

global s  % tension parameter
s = 0.5;

global nsurrogates; % number of surrogates for bootstrapping
nsurrogates = 10000;


%model_coefficients = 0.07*[hann(20)', -0.5*ones(20,1)'];


%model_coefficients = 0.07*hann(20)';
% model_coefficients = [1.8507
%     -0.6429
%     -0.4042
%     0.1227
%     0.0657
%     -0.0658
%     0.0196
%     0.1035
%     0.0058
%     -0.0971
%     -0.0445
%     0.0307
%     0.0280
%     0.0097
%     0.0097
%     0.0130
%     7.9320e-04
%     -0.0202
%     -0.0230
%     0.0192]';
model_coefficients = [-0.2314 %%%%% spectral peak, high freq
    0.1613
    0.1405
    0.0741
    -0.0098
    -0.0836
    -0.1193
    -0.1048
    -0.0585
    0.0011
    0.0559
    0.0874
    0.0837
    0.0614
    0.0289
    -0.0050
    -0.0316
    -0.0423
    -0.0285
    0.0185]';




nlags = size(model_coefficients,2);

nobs = 1000;   % number of observations per trial
sgnl = 1;
nelectrodes = 1;

f0 = 500;  % sampling frequency (Hz)
T = nobs/f0;      % total length of recording (seconds)
dt = 1/f0; % seconds

N = T*f0 +nlags;   % number of samples needed +2
df = 1/T;   % frequency resolution
fNQ = f0/2; % Nyquist frequency

taxis = dt:dt:T; % time axis
noise = 0.25;

model_order =30; % order used in model estimation
adj_true = 1;



%%%  Generate data ---------------------------------------------------

b = zeros(1,1,nlags);
b(1,1,:)= model_coefficients;


%%% Generate white noise data from above coefficients
data = zeros(1,N);
i=1;
for k = nlags:length(data)-1
    data(:,k+1) = myPrediction(data(:,1:k),b);
    noise_process(i) = noise.*randn(size(data,1),1);
    data(:,k+1) = data(:,k+1) + noise_process(i);
    i=i+1;
end

data = data(:,nlags+1:end);



%%%  Plot data ---------------------------------------------------



% subplot(2,2,[1 2])
% plot(taxis,data(1,:));
% 
% ylabel('Signal')
% xlabel('Time (seconds)')
% 
% title('Simulated Signal','FontSize',15);

ax1 = subplot(2,2,3);
mySpec(data(1,:),f0,'yesplot','tapers');
title('True Signal Spectrogram','FontSize',15);

%%% Fit spline to data ---------------------------------------------------
%cntrl_pts = [0:1:model_order]
cntrl_pts = make_knots(model_order,floor(model_order/3));
adj_standard = build_ar(data,model_order);
[ adj_mat] = build_ar_splines( data, cntrl_pts(end), cntrl_pts );
[bhat, yestimate] = estimate_coefficient_fits( data, adj_mat,  cntrl_pts(end),cntrl_pts );
[b_est_stand, y_est_stand] = estimate_standard( data, adj_standard, model_order );
%model_order = cntrl_pts(end);
%%% Plot results ---------------------------------------------------------

ax2=subplot(2,2,4);
mySpec(yestimate,f0,'yesplot','tapers');
linkaxes([ax1,ax2],'y')
title('Estimated Signal Spectrogram','FontSize',15);

subplot(2,2,[1,2])
    plot(taxis(model_order+1:end),data(:,model_order+1:end),'k','LineWidth',1.5);
    hold on;
    plot(taxis(model_order+1:end),yestimate,'--r','LineWidth',2);
   title('Single Node Simulation','FontSize',20)
xlabel('Time (seconds)','FontSize',17)
h = legend('True Signal','Spline Estimated Signal');
set(h,'FontSize',15);

%%% Plot all results --------------------------------------------
goodness_of_fit_residuals;
goodness_of_fit_bootstrap;
goodness_of_fit_spectrum;
%%%% Make table of all inputs ------------------------
% Sampling_Frequency = f0;
% Noise_Variance = noise.^2;
% T_seconds = T;
% True_Model_Order = nlags;
% Estimated_Order = model_order;
% DW_statistic = dw;
% Tension_Parameter = s;
% Tp = table(Sampling_Frequency,Noise_Variance,T_seconds,True_Model_Order,...
%     Estimated_Order,Tension_Parameter,moAIC,moBIC, DW_statistic, pval);
% figure;
% uitable('Data',Tp{:,:}','RowName',Tp.Properties.VariableNames,...
%     'Units', 'Normalized', 'Position',[0,0,1,1]);


% % Save all simulation and table plots ---------------------------

% 
% h = get(0,'children');
% for i=1:length(h)
%     saveas(h(i), ['1N_ReverseEngineered_modelorder_' num2str(model_order) '_' num2str(i)], 'jpg');
% end