%%% ar_real_data
clear all
close all
f0 = 499.9;
load('data_unfiltered.mat')
%load('data.mat')
data = data - repmat(mean(data,2), [1,129290]);
data_type = 'real';

N = 1000;
sz = 9e4:9e4+N-1;
pre_sz = 2e4:2e4+N-1;
ntwk=34;
% badchannels = [1,9,21,32,83, 8,31];
ntwk = [2 4 6 18 22 34 36 38 57 64];
ntwk = [2 6 18 22 42 46 90 82 77 ];
data_sz= data(ntwk,sz);
data_presz= data(ntwk,pre_sz);

%y=fft(data(9,:));

% X = data_presz;
% mvar_aic;  
% % moAIC = 6; moBIC =3 (10 seconds) %%% for 10 node network
% % moAIC = 3; moBIC =2 (2.5 seconds)
% 
% X = data_sz;
% mvar_aic;
% moAIC = 10; moBIC =4 (10 seconds)
% moAIC = 4; moBIC =3 (2.5 seconds)

%%
figure 
subplot 211
plotchannels(data_presz')
title('Channels Pre Seizure')
subplot 212
plotchannels(data_sz')
title('Channels During Seizure')
% h = get(0,'children');
% saveas(h, ['realdata_channels'], 'jpg');
% close all

% figure;
% load('ECoG.mat')
% plot(ECoG.ElectrodeXY(:,1),ECoG.ElectrodeXY(:,2),'.')
% text(ECoG.ElectrodeXY(:,1),ECoG.ElectrodeXY(:,2),num2str((1:96)'))

%%%% ---- Define inputs
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

nlags = 20;
ctpts = make_knots(nlags,floor(nlags/3)); %0:5:nlags;
%%
%%%% --------------- PRE SEIZURE ------------------------------------
% X = data_presz;
% mvar_aic; 
morder=20;
tic
[ adj_standard ] = build_ar( data_presz, morder); % Build network using splines
standard_time  = toc;
fprintf('\nMVGC time = %d\n',standard_time);
fprintf('using model order = %d\n',morder);
tic
[ adj_spline ] = build_ar_splines( data_presz, nlags,ctpts); % Build network using splines
spline_time  = toc;
fprintf('\nspline MVGC time = %d\n',spline_time);
fprintf('using model order = %d\n',nlags);
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

suptitle('Pre-Seizure')



% h = get(0,'children');
% saveas(h, ['realdata_preseizure_network'], 'jpg');
% close all


[ bhat, yestimate ] = estimate_coefficient_fits( data_presz, adj_spline, nlags, ctpts);
model_order = nlags;
nelectrodes = length(ntwk);
adj_mat = adj_spline;
cntrl_pts = ctpts;
data = data_presz;
adj_true = adj_spline;
b= zeros(nelectrodes,nelectrodes,model_order);
%goodness_of_fit_residuals;
%goodness_of_fit_bootstrap;

% h = get(0,'children');
% j=1;
% for i=length(h):-1:1
%     saveas(h(j), ['realdata_preseizure_plots'   num2str(i)], 'jpg');
%     j=j+1;
% end
% close all


%%
%%%% --------------- SEIZURE ------------------------------------
% X = data_sz;
% mvar_aic; 
morder=50;
% tic
% [ adj_standard ] = build_ar( data_sz, morder); % Build network using splines
% standard_time  = toc;
% fprintf('\nMVGC time = %d\n',standard_time);
% fprintf('using model order = %d\n',morder);
nlags = 50; 
% ctpts = [0 5 10 15 20 23 25 26 27 28 29 30 31 32 33 34 35 36 38 40 42 45 47 49 50];
ctpts =   [0 2 4 10 15 20 25 30 33 35 35 40 50];
tic
[ adj_spline ] = build_ar_splines( data_sz, nlags,ctpts); % Build network using splines
spline_time  = toc;
fprintf('\nspline MVGC time = %d\n',spline_time);
fprintf('using model order = %d\n',nlags);

figure();
subplot 221
%plotNetwork(adj_standard) 
%title(strcat({'Standard, '},num2str(standard_time),{' s'}))
% subplot 222
% plotGraph(adj_standard)
subplot 223
plotNetwork(adj_spline)  % plot network
title(strcat({'Spline, '},num2str(spline_time),{' s'}))
colorbar
subplot 224
plotGraph(adj_spline)

suptitle('Seizure')

% h = get(0,'children');
% saveas(h, ['realdata_seizure_network'], 'jpg');
% close all

[ bhat, yestimate ] = estimate_coefficient_fits( data_sz, adj_spline, nlags, ctpts);
model_order = nlags;
nelectrodes = length(ntwk);
adj_mat = adj_spline;
cntrl_pts = ctpts;
data = data_sz;
adj_true = adj_spline;
b= zeros(nelectrodes,nelectrodes,model_order);
% 
% goodness_of_fit_residuals;
 goodness_of_fit_bootstrap;
% 
% 
% h = get(0,'children');
% j=1;
% for i=length(h):-1:1
%     saveas(h(j), ['realdata_seizure_plots'   num2str(i)], 'jpg');
%     j=j+1;
% end
% close all



