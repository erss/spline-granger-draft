%%%%%%%% simulate data --------------------------------------------------

T = 5;      % total length of recording (seconds)
dt = 0.001; % seconds

f0 = 1/dt;  % sampling frequency (Hz)
N1 = T*f0;   % number of samples needed
df = 1/T;   % frequency resolution
fNQ = f0/2; % Nyquist frequency

N = T*f0 + 40;
taxis = dt:dt:T; % time axis
noise = 7;
data = zeros(3,N);


%%% SIM 1
a1 = 0.07*[hann(20)', -0.5*ones(20,1)']';   % AR coefficients for signal 1
a2 = 0.05*[-0.5*ones(20,1)', hann(20)']';   %                  ...signal 2
a3 = -.3*ones(size(a1));                    %                  ...signal 3
a = zeros(3,3,40);
a(1,1,:) = a1;
a(1,2,:) = a2;
a(2,2,:) = a2;
a(3,3,:) = a3;
% a(:,:,1) = [a1';a2';zeros(size(a2'))];
% a(:,:,2) = [zeros(size(a2'));a2';zeros(size(a2'))];
% a(:,:,3) = [zeros(size(a2'));zeros(size(a2'));a3'];

nlags = length(a1);


for k = nlags:length(data)-1;
    data(:,k+1) = myPrediction(data(:,1:k),a);
    data(:,k+1) = data(:,k+1) + noise.*randn(size(data,1),1);
end

mvar_aic;



subplot(3,2,[1 2])
 plot(data(1,:));
 hold on;
 plot(data(2,:));
 plot(data(3,:));

ylabel('Signal')
xlabel('Time (seconds)')
legend('x1','x2','x3')
title('Simulated Signal','FontSize',15);


subplot(3,2,3)
mySpec(data(1,:),f0);

%%% Fit spline to data ---------------------------------------------------
model_order =40;
flag = 1; % use splines
cntrl_pts = make_knots(model_order,10);
%cntrl_pts = 0:10:nlags;
[ adj_mat] = build_ar_splines( data, model_order, cntrl_pts );
%%
%[ bhat, yhat ] = estimate_coef( data, adj_mat, nlags, flag,cntrl_pts);
[ bhat, yhat ] = estimate_coefficient_fits( data, adj_mat, nlags, cntrl_pts);

%%


%%
subplot(3,2,[5 6])
% plot(bhat,'LineWidth',1.5)
% hold on
% plot(cntrl_pts(2:end),bhat(cntrl_pts(2:end)),'o')
% title('Estimated Coefficients','FontSize',15);
% figure;
% plot(data);
% hold on
% plot(yhat,'--r');


subplot(3,2,4)
mySpec(yhat(1,:),f0);
title('Estimated signal spectrogram','FontSize',15);

figure; plotSignals(data)
figure; plotSignals(yhat)

 figure; 
for i = 1:3
    for k = 1:3
      
plot(squeeze(bhat(i,k,:)),'--r','LineWidth',1.5);
hold on;
plot(cntrl_pts(2:end),squeeze(bhat(i,k,cntrl_pts(2:end))),'ro')
plot(squeeze(a(i,k,:)),'k','LineWidth',1.5);
    end

end

