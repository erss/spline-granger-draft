%%% 1 Node network simulations, AR(2)
%% AR(1) Low Frequency
T = 5;      % total length of recording (seconds)
dt = 0.001; % seconds

f0 = 1/dt;  % sampling frequency (Hz)
N1 = T*f0;   % number of samples needed
df = 1/T;   % frequency resolution
fNQ = f0/2; % Nyquist frequency

taxis = dt:dt:T; % time axis


%a1 = [0.9];   % AR(1) Low Frequency
%a1 = [0.9; -0.8];   % AR(2) high freq
%a1 = [0.3; 0.3]; % AR(2) Low Frequency
%a1 = [0.9; -0.1];  % AR(2) Low frequency
a1 = [0.5];
L = length(a1);                             % Number of AR terms.
N = N1+L;                                 % Number of time steps.

x1 = zeros(N,1);    

noise = 7;                                 % noise level

for k=L+1:N                                % For each time step,
    x1(k) = sum(a1.*x1(k-1:-1:k-L)) + noise*randn();
   %...x(now) depends on past.
end
x1 = x1(L+1:end);                          % Drop the first L terms (0s).   


figure;
subplot (2,3,[1 3])                              %Plot signals.
plot(taxis, x1)

ylabel('Signal','FontSize',15)
xlabel('Time (seconds)','FontSize',15)
title(num2str(a1'),'FontSize',15);

subplot(2,3,[4 5])
plot(taxis(1:125),x1(1:125))
axis tight
ylabel('Signal','FontSize',15)
xlabel('Time (seconds)','FontSize',15)
title('Zoomed signal','FontSize',15);

subplot(2,3,6)
mySpec(x1,f0);

%%% Define inputs for model building and GoF -------------------------
data = x1';                      % Conglomerate data in one matrix
nlags = 200;                               % Define order of AR model
                                           % needs to be larger than true
                                           % order

true_coeffs = cell(1);
true_coeffs{1,1} = padarray(a1,nlags-L,'post'); % True coefficents up to 
                                                %     number of lags being estimated

tic
[ adj_spline ] = build_ar_splines( data, nlags); % Build network using splines
splinetime  = toc;
tic
[ adj_stand] = build_ar( data, nlags );          % Build network using standard ar
standardtime = toc;

which_electrode = 1;
b_est = ar_gof(adj_spline, adj_stand,data, nlags,true_coeffs,taxis,which_electrode); % Run goodness of fit tests




%% Low frequency, bandpass filtered

T = 5;      % total length of recording (seconds)
dt = 0.001; % seconds

f0 = 1/dt;  % sampling frequency (Hz)
N1 = T*f0;   % number of samples needed
df = 1/T;   % frequency resolution
fNQ = f0/2; % Nyquist frequency

taxis = dt:dt:T; % time axis


a1 = [0.3];   %AR coefficients for signal 1

L = length(a1);                             % Number of AR terms.
N = N1+L;                                 % Number of time steps.

x1 = zeros(N,1);    

noise = 1;                                 % noise level

for k=L+1:N                                % For each time step,
    x1(k) = sum(a1.*x1(k-1:-1:k-L)) + noise*randn();
   %...x(now) depends on past.
end
x1 = x1(L+1:end);                          % Drop the first L terms (0s).   


figure;
subplot (2,3,[1 3])                              %Plot signals.
plot(taxis, x1)

ylabel('Signal','FontSize',15)
xlabel('Time (seconds)','FontSize',15)
title('AR(1)','FontSize',15);

subplot(2,3,[4 5])
plot(taxis(1:125),x1(1:125))
axis tight
ylabel('Signal','FontSize',15)
xlabel('Time (seconds)','FontSize',15)
title('Zoomed signal','FontSize',15);

subplot(2,3,6)
mySpec(x1,f0);

order = 300;
band = [10 20]; % if BAND(1) >0 and BAND(2) <fNQ
f = [0 (band(1)-1)/fNQ band/fNQ (band(2)+1)/fNQ 1];
a = [0 0 1 1 0 0];


flt = firls(order,f,a);
y = filtfilt(flt,1,x1);

figure;

subplot (2,3,[1 3])                              %Plot signals.
plot(taxis, y)

ylabel('Signal','FontSize',15)
xlabel('Time (seconds)','FontSize',15)
title('Bandpass filtered Low Freq','FontSize',15);

subplot(2,3,[4 5])
plot(taxis(1:125),y(1:125))
axis tight
ylabel('Signal','FontSize',15)
xlabel('Time (seconds)','FontSize',15)
title('Zoomed signal','FontSize',15);

subplot(2,3,6)
mySpec(y,f0);





%%% Define inputs for model building and GoF -------------------------
data = y';                      % Conglomerate data in one matrix
nlags = 200;                               % Define order of AR model
                                           % needs to be larger than true
                                           % order

true_coeffs = cell(1);
true_coeffs{1,1} = padarray(a1,nlags-L,'post'); % True coefficents up to 
                                                %     number of lags being estimated

tic
[ adj_spline ] = build_ar_splines( data, nlags); % Build network using splines
splinetime  = toc;
tic
[ adj_stand] = build_ar( data, nlags );          % Build network using standard ar
standardtime = toc;
% 
% subplot 234; 
% adj_true = [1];          % Plot true network
% plotNetwork(adj_true)
% title('Truth');
% subplot 235;                               % Plot network using standard regressors
% str = {num2str(standardtime), ' seconds' };
% plotNetwork(adj_stand)
% title(strcat({'Standard, '},num2str(standardtime),{' s'}))
% subplot 236;                               % Plot network using spline regressors
% plotNetwork(adj_spline)
% title(strcat({'Spline, '},num2str(splinetime),{' s'}))

which_electrode = 1;
b_est = ar_gof(adj_spline,adj_stand, data, nlags,true_coeffs,taxis,which_electrode);
 
