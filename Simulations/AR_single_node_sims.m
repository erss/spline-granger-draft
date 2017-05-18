%%%%%%%% Single node network simulations ----------------------------------
clear all;

noise_type = 'pink';   % 'white' or 'pink'
frequency_type = 'low'; % 'low' or 'high'


nlags = 0;
n=3;
if strcmp(noise_type,'white')
    n=2;
    if strcmp(frequency_type,'high') % ----high frequency
        model_coefficients = [0.9 -0.8];
        nlags = size(model_coefficients,2);
    elseif strcmp(frequency_type,'low') % ----low frequency
         % model_coefficients = 0.9;
          model_coefficients = [0.3 0.3];
         % model_coefficients = [0.9 -0.1];
         nlags = size(model_coefficients,2);
    end
end


%%% Define model inputs ---------------------------------------------------
nobs = 1000;   % number of observations per trial
sgnl = 1;
nelectrodes = 1;

f0 = 200;  % sampling frequency (Hz)
T = nobs/f0;      % total length of recording (seconds)
dt = 1/f0; % seconds

N = T*f0 +nlags;   % number of samples needed +2
df = 1/T;   % frequency resolution
fNQ = f0/2; % Nyquist frequency

taxis = dt:dt:T; % time axis
noise = 0.25;
data = zeros(1,N);

model_order = 20; % order used in model estimation

%%%  Generate data ---------------------------------------------------
if strcmp(noise_type,'white') % ------------ WHITE NOISE -------

    
    b = zeros(1,1,nlags);
    b(1,1,:)= model_coefficients;
      

    %%% Generate white noise data from above coefficients
    i=1;
    for k = nlags:length(data)-1;
       data(:,k+1) = myPrediction(data(:,1:k),b);
       noise_process(i) = noise.*randn(size(data,1),1);
       data(:,k+1) = data(:,k+1) + noise_process(i);
          i=i+1;
    end
    
    data = data(:,nlags+1:end);
    [faxis,S] = myTheoreticalSpectrum(model_coefficients,noise_process,f0);
    

elseif strcmp(noise_type,'pink') % ------------ PINK NOISE -------
 alpha = 0.33;
 data  = make_pink_noise(alpha,nobs,dt);
end

%%%  Plot data ---------------------------------------------------



subplot(n,2,[1 2])
 plot(taxis,data(1,:));

ylabel('Signal')
xlabel('Time (seconds)')
legend('x1')
title('Simulated Signal','FontSize',15);

subplot(n,2,3)
mySpec(data(1,:),f0);

if strcmp(noise_type,'white')
    hold on
    plot(faxis,S,'r','LineWidth',1.5);
end
%%% Fit spline to data ---------------------------------------------------

cntrl_pts = make_knots(model_order,10);
[ adj_mat] = build_ar_splines( data, model_order, cntrl_pts );
[bhat, yhat] = estimate_coefficient_fits( data, adj_mat, model_order,cntrl_pts );

%%% Plot results ---------------------------------------------------------



if strcmp(noise_type,'pink')
    subplot(n,2,[5 6])
    plot(dt:dt:(model_order/f0),squeeze(bhat(1,1,:)),'LineWidth',1.5)
    hold on
    plot(cntrl_pts(2:end)./f0,squeeze(bhat(1,1,cntrl_pts(2:end))),'o')
    title('Estimated Coefficients','FontSize',15);

end

subplot(n,2,4)
mySpec(yhat,f0);

if strcmp(noise_type,'white')
    hold on
    plot(faxis,S,'r','LineWidth',1.5);
end


title('Estimated signal spectrogram','FontSize',15);


%%% Determine what AIC thinks is best order

if strcmp(noise_type,'white')
   % mvar_aic;
    goodness_of_fit_bootstrap;
end

goodness_of_fit_spectrum;



