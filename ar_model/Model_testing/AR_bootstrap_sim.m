%%% Generate data -------------------------------------------------------
T = 5;      % total length of recording (seconds)
dt = 0.001; % seconds
nelectrodes = 3;
f0 = 1/dt;  % sampling frequency (Hz)
N1 = T*f0;   % number of samples needed
df = 1/T;   % frequency resolution
fNQ = f0/2; % Nyquist frequency

taxis = dt:dt:T; % time axis


a1 = 0.07*[hann(20)', -0.5*ones(20,1)']';   %AR coefficients for signal 1
a2 = 0.05*[-0.5*ones(20,1)', hann(20)']';    %                  ...signal 2
a3 = -.3*ones(size(a1));                    %                  ...signal 3

L = length(a1);                             % Number of AR terms.
N = N1+L;                                 % Number of time steps.

x1 = zeros(N,1);
x2 = zeros(N,1);                           %Create the AR model,
x3 = zeros(N,1);    

noise = 7;                                 % noise level

for k=L+1:N                                % For each time step,
    x1(k) = sum(a1.*x1(k-1:-1:k-L)) + sum(a2.*x2(k-1:-1:k-L)) +  noise*randn();
    x2(k) = sum(a2.*x2(k-1:-1:k-L)) + noise*randn();
    x3(k) = sum(a3.*x3(k-1:-1:k-L)) + noise*randn();  %...x(now) depends on past.
end
x1 = x1(L+1:end);                          % Drop the first L terms (0s).   
x2 = x2(L+1:end);                             
x3 = x3(L+1:end);


%%% Define inputs for model building -----------------------------------
data = [x1';x2';x3'];                      % Conglomerate data in one matrix
nlags = 40;                               % Define order of AR model
                                           % needs to be larger than true
                                           % order
model_order = nlags;
%%% Build network using splines -----------------------------------------
cntrl_pts = make_knots(nlags,10);
[ adj_spline ] = build_ar_splines( data, nlags,cntrl_pts); 
b = estimate_coef(data,adj_spline, nlags,1); % estimated coefficients  from adj_spline

%%% Bootstrap for each coefficient estimate on the surrogate data ---------

electrode = 1; % which electrode in network to generate data for 
nsurrogates = 1000; % number of surrogates

[b_surrogates, bounds]= myBootstrap(data,adj_spline,nlags,electrode,nsurrogates,cntrl_pts);

   j=1;
      for p = 1:nelectrodes
        if adj_mat(electrode,p) == 1
            bhat(electrode,p,:) = b(j:j+model_order-1);
            j= j+model_order;
        end    
     end



figure;
for row = 1:nsurrogates
 plot(b_surrogates(row,:))  % plot surrogate data
 hold on;
end
figure;
plot(b(electrode,:),'--k','LineWidth',1) % plot original model fit


%%% plot estimated coefficients and confidence bounds --------------------
figure;
plot(b(electrode,:),'r');
hold on
plot((bounds(1,:)),'--r');
plot((bounds(2,:)),'--r');
plot([a1' a2' zeros(1,40)],'k')



