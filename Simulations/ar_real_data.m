%%% spectral confidence bounds, real data
%%% Spectral confidence bounds

% load('MG49_Seizure36.mat', 'ECoG') 
% load('ECoG.mat')
% 
% % Reference data
% 
% 
% re_referenced_data = bipolarRefSimulation(ECoG.Data);                        %Re-reference the data,
% 
% BAND = [4  50]; 
% fs = ECoG.SamplingRate;                               % Save the sampling rate.
% filtered_data = lsfilter(re_referenced_data, fs, BAND);        % Filter the data,
%               
% %data = re_referenced_data';
% data = filtered_data';
% save('data');
%%
clear all
load('data_unfiltered.mat')
data_type = 'real';
sz = 9e4:9e4+5000-1;
pre_sz = 2e4:2e4+4999;

%ntwk = [2 14 48 60 34];
ntwk = [2 4 6 18 22 34 36 38 57 64];

data_sz= data(ntwk,sz);
data_presz= data(ntwk,pre_sz);

figure 
subplot 211
plotchannels(data_presz')
title('Channels Pre Seizure')
subplot 212
plotchannels(data_sz')
title('Channels During Seizure')

% figure;
% load('ECoG.mat')
% plot(ECoG.ElectrodeXY(:,1),ECoG.ElectrodeXY(:,2),'.')
% text(ECoG.ElectrodeXY(:,1),ECoG.ElectrodeXY(:,2),num2str((1:96)'))
%%
N = 5000;%129290;    % number of samples 
f0 = 499.7;     % sampling frequency (Hz)

T = N/f0;         % total length of recording (seconds)
dt = 1/f0;    % seconds

df = 1/T;      % frequency resolution
fNQ = f0/2;    % Nyquist frequency
noise= 0.25;
taxis = dt:dt:T; % time axis

global s;
s = 0.5;
global nsurrogates;
nsurrogates = 100;

nlags = 100;
ctpts = make_knots(100,10);
tic
%[ adj_standard ] = build_ar( data_presz, nlags); % Build network using splines
standard_time  = toc;
tic
[ adj_spline ] = build_ar_splines( data_presz, nlags,ctpts); % Build network using splines
spline_time  = toc;

figure();
subplot 221
%plotNetwork(adj_standard) 
title(strcat({'Standard, '},num2str(standard_time),{' s'}))
colorbar
subplot 222
%plotGraph(adj_standard)
subplot 223
plotNetwork(adj_spline)  % plot network
title(strcat({'Spline, '},num2str(spline_time),{' s'}))
colorbar
subplot 224
plotGraph(adj_spline)

suptitle('Pre-Seizure')

[ bhat, yestimate ] = estimate_coefficient_fits( data_presz, adj_spline, nlags, ctpts);
model_order = nlags;
nelectrodes = length(ntwk);
adj_mat = adj_spline;
cntrl_pts = ctpts;
data = data_presz;
adj_true = adj_spline
b= zeros(nelectrodes,nelectrodes,model_order);
goodness_of_fit_residuals;
%goodness_of_fit_bootstrap;
%goodness_of_fit_spectrum


%%

N = 5000;%129290;    % number of samples 
f0 = 499.7;     % sampling frequency (Hz)

T = N/f0;         % total length of recording (seconds)
dt = 1/f0;    % seconds

df = 1/T;      % frequency resolution
fNQ = f0/2;    % Nyquist frequency
noise= 0.25;
taxis = dt:dt:T; % time axis

global s;
s = 0.5;
global nsurrogates;
nsurrogates = 100;

nlags = 100;
ctpts = make_knots(100,10);
tic
[ adj_standard ] = build_ar( data_sz, nlags); % Build network using splines
standard_time  = toc;
tic
[ adj_spline ] = build_ar_splines( data_sz, nlags,ctpts); % Build network using splines
spline_time  = toc;

figure();
subplot 221
plotNetwork(adj_standard) 
title(strcat({'Standard, '},num2str(standard_time),{' s'}))
subplot 222
plotGraph(adj_standard)
subplot 223
plotNetwork(adj_spline)  % plot network
title(strcat({'Spline, '},num2str(spline_time),{' s'}))
colorbar
subplot 224
plotGraph(adj_spline)

suptitle('Seizure')
%%
[ bhat, yhat ] = estimate_coefficient_fits( data_sz, adj_spline, nlags, ctpts);
model_order = nlags;
nelectrodes = length(ntwk);
adj_mat = adj_spline;
cntrl_pts = ctpts;
data = data_sz;
adj_true = adj_spline;
b= zeros(nelectrodes,nelectrodes,model_order);
goodness_of_fit_bootstrap;
%goodness_of_fit_spectrum

%%
% 
% %%
% nelectrodes = size(data,1);
% electrode = randi(nelectrodes);
% nrealizations = 1;           
% 
% 
% %y = data(electrode,:);   
% figure(2);
% subplot 131
% [faxis, h] = mySpec( data(electrode,:), f0 );   % plot spectrum of true
% title('ECoG Spectrum');
%    
% 
%  
%  flag = 1; % use splines to estimate
% [bhat, yhat] = estimate_coef(data,adj_spline,nlags,flag);
% figure(8);
% 
% residuals_spline_model = data(electrode,nlags+1:end)-yhat(electrode,:);
% 
% 
% subplot 131
% plot(residuals_spline_model);
% title('residuals spline model');
% 
% 
% subplot 132
% autocorr(residuals_spline_model);
% 
% 
% subplot 133
% n= length(residuals_spline_model);
% [F,x] = ecdf(residuals_spline_model);
% Fn = normcdf(x,mean(residuals_spline_model),std(residuals_spline_model));
% 
% plot(F,F,'b','LineWidth',2);
% hold on;
% plot(F,F-1.36/sqrt(n),'--r','LineWidth',2);
% plot(F,F+1.36/sqrt(n),'--r','LineWidth',2);
% plot(F,Fn,'g','LineWidth',2);
% 
% legend('1-1 line', 'Lower CI', 'Upper CI','Emp vs Theoretical');
% xlim([0 1]);
% ylim([0 1]);
% 
% for k =1:size(yhat,1)
%     figure(3);
%     subplot (nelectrodes/5,5,k)
%    plot(yhat(k,:),'r')   % plot all glmfit outputs
%    figure(4);
%    subplot (nelectrodes/5,5,k)
%    plot(data(k,:),'k')  % plot ECoG signals
% end
% figure(3)
% suptitle('glmfit signals')
% figure(4)
% suptitle('ECoG signals')
% 
% 
% 
% h_sum = 0;
% for i = 1:nrealizations
%         data_hat = zeros(nelectrodes,N);
%     for k = nlags:length(data_hat)-1;
%         data_hat(:,k+1) = myPrediction(data_hat(:,1:k),bhat,nlags);
%         data_hat(:,k+1) = data_hat(:,k+1) + noise.*randn(nelectrodes,1);
%     end
%     data_hat= data_hat(:,41:end);
%     y0 = data_hat(electrode,:);   
%     [faxis, h_hat] = mySpec( y0, f0,0 ); % compute spectra
%     h_sum = h_hat + h_sum;
% end
% h_hat = h_sum/nrealizations;
% 
% 
% figure(5);
% for k =1:size(data_hat,1)
%     subplot (nelectrodes/5,5,k)
%    plot(taxis(41:end),data_hat(k,:)) 
% end
% suptitle('Simulated data using bhat');
% 
% 
% figure(2)
% subplot 132
% 
% plot(faxis,h_hat);     
% xlim([0 f0/4]);
% xlabel('Frequency (Hz)','FontSize',15);
% ylabel('Power','FontSize',15);
% title('simulated signal','FontSize',15);
% 
% subplot 133
% [faxis, hglm] = mySpec(yhat(electrode,:),f0);
% title('glmfit signal');
% 
% % compute and plot cumulative distributions of spectra
% [H, X] = ecdf(h);           
% [H1, X1] = ecdf(h_hat);
% [H2, X2] = ecdf(hglm);
% 
% % 
% % H=2*H;
% % H1 = 2*H1;
% % H2 = 2*H2;
% 
% 
% figure(6);
% plot(X2,H2,'r','LineWidth',1.5);
% hold on
% plot(X,H,'k','LineWidth',1.5);
% legend('Estimated Signal','True Signal')
% title('CDFs of Spectrum');
% 
% 
% 
% % Compute confidence bounds for estimated signal (Priestley p 478)
% 
% a = 2.2414; % for 95% confidence bounds
% N = size(data,2); % number of observations from which H is computed ? length of signal ?
% 
% flag = 'biased'; % divide by 1/N
% 
% R = xcov(yhat(electrode,:),flag); % autocovariance of estimated signal
% 
% G = sum(R(3:end-2).^2);
% G = G/(4*pi);
% 
% conf1 = a*sqrt(8*pi*G/N);
% 
% %plot(X2,H2 + conf1, '--r');
% %plot(X2,H2 - conf1, '--r');
% 
% 
% % figure;
% % myKS(data(electrode,:),yhat(electrode,:))
% % figure;
% % myKS(data(electrode,:),data_hat(electrode,:))
