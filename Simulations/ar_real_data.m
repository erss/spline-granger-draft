%%% ar_real_data
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
h = get(0,'children');
saveas(h, ['realdata_channels'], 'jpg');
close all

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


%%%% --------------- PRE SEIZURE ------------------------------------

nlags = 100;
ctpts = make_knots(100,10);
tic
[ adj_standard ] = build_ar( data_presz, nlags); % Build network using splines
standard_time  = toc;
tic
[ adj_spline ] = build_ar_splines( data_presz, nlags,ctpts); % Build network using splines
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

suptitle('Pre-Seizure')


h = get(0,'children');
saveas(h, ['realdata_preseizure_network'], 'jpg');
close all


[ bhat, yestimate ] = estimate_coefficient_fits( data_presz, adj_spline, nlags, ctpts);
model_order = nlags;
nelectrodes = length(ntwk);
adj_mat = adj_spline;
cntrl_pts = ctpts;
data = data_presz;
adj_true = adj_spline
b= zeros(nelectrodes,nelectrodes,model_order);
goodness_of_fit_residuals;
goodness_of_fit_bootstrap;

h = get(0,'children');
j=1;
for i=length(h):-1:1
    saveas(h(j), ['realdata_preseizure_plots'   num2str(i)], 'jpg');
    j=j+1;
end
close all



%%%% --------------- SEIZURE ------------------------------------
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

h = get(0,'children');
saveas(h, ['realdata_seizure_network'], 'jpg');
close all

[ bhat, yestimate ] = estimate_coefficient_fits( data_sz, adj_spline, nlags, ctpts);
model_order = nlags;
nelectrodes = length(ntwk);
adj_mat = adj_spline;
cntrl_pts = ctpts;
data = data_sz;
adj_true = adj_spline;
b= zeros(nelectrodes,nelectrodes,model_order);

goodness_of_fit_residuals;
goodness_of_fit_bootstrap;


h = get(0,'children');
j=1;
for i=length(h):-1:1
    saveas(h(j), ['realdata_seizure_plots'   num2str(i)], 'jpg');
    j=j+1;
end
close all



