%%%%%%%% Three node network simulations -----------------------------------
clear all;
%%% Define model inputs ---------------------------------------------------

nelectrodes = 3; % number of electrodes
nlags = 40;  % true model order
model_order = 100; % order used in model estimation
T = 5;      % total length of recording (seconds)
dt = 0.001; % seconds

f0 = 1/dt;  % sampling frequency (Hz)
df = 1/T;   % frequency resolution
fNQ = f0/2; % Nyquist frequency

N = T*f0 + nlags;
taxis = dt:dt:T; % time axis
noise = 0.25;
data = zeros(3,N);


%%% SIM 1 ----------------------------------------------------------------
a1 = 0.07*[hann(20)', -0.5*ones(20,1)']';   
a2 = 0.05*[-0.5*ones(20,1)', hann(20)']';   
a3 = -.3*ones(size(a1));                   

a = zeros(3,3,40);                 % Model coefficients                 
a(1,1,:) = a1;            
a(1,2,:) = a2;                                               
a(2,2,:) = a2;
a(3,3,:) = a3;
nlags = length(a1);
adj_true = [1 1 0; 0 1 0; 0 0 1];  % True network stucture

%%% SIM 2 ----------------------------------------------------------------
% a1 = 0.07*[hann(20)', -0.5*ones(20,1)']';   
% a2 = 0.03*[-0.5*ones(20,1)', hann(20)']';   
% a3 = -.3*ones(size(a1));
% a = zeros(3,3,40);                         % Model coefficients
% a(1,1,:) = a1;
% a(1,2,:) = a2;
% a(2,1,:) = a1;
% a(2,2,:) = a2;
% a(3,3,:) = a3;
% adj_true = [1 1 0; 1 1 0; 0 0 1];          % True network structure

%%% SIM 3----------------------------------------------------------------
% a1 = 0.07*[hann(20)', -0.5*ones(20,1)']';   
% a2 = 0.03*[-0.5*ones(20,1)', hann(20)']';
% a = zeros(3,3,40);
% a(1,3,:) = a1;                             % Model coefficients
% a(2,3,:) = a1;
% a(3,3,:) = a2;
% adj_true = [0 0 1; 0 0 1; 0 0 1];          % True network structure




%%% Simulate data ------------------------------------------
for k = nlags:length(data)-1;
    data(:,k+1) = myPrediction(data(:,1:k),a);
    data(:,k+1) = data(:,k+1) + noise.*randn(size(data,1),1);
end

% mvar_aic; run to see mvgc toolbox order result



subplot(2,3,[1 3])
 plot(data(1,:));
 hold on;
 plot(data(2,:));
 plot(data(3,:));

ylabel('Signal')
xlabel('Time (seconds)')
legend('x1','x2','x3')
title('Simulated Signal','FontSize',15);

% 
% subplot(2,2,3)
% mySpec(data(1,:),f0);

%%% Fit standard AR to data ----------------------------------------------

tic
[ adj_standard] = build_ar( data, model_order);
standardtime  = toc;
%%% Fit spline to data ---------------------------------------------------

cntrl_pts = make_knots(model_order,10);
tic
[ adj_mat] = build_ar_splines( data, model_order, cntrl_pts );
splinetime  = toc;
[ bhat, yhat ] = estimate_coefficient_fits( data, adj_mat, model_order, cntrl_pts);


%%% Plot results ----------------------------------------------------------
% 
% subplot(2,2,4)
% mySpec(yhat(1,:),f0);
% title('Estimated signal spectrogram','FontSize',15);


subplot(2,3,4)
plotNetwork(adj_true)
title('True Network')

subplot(2,3,5)
plotNetwork(adj_standard)
title(strcat({'Standard, '},num2str(standardtime),{' s'}))


subplot(2,3,6)
plotNetwork(adj_mat)
title('Spline Network')
title(strcat({'Spline, '},num2str(splinetime),{' s'}))


b=a;
goodness_of_fit_spectrum;
goodness_of_fit_bootstrap;