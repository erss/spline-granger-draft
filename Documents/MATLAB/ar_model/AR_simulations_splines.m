%% Simulation 1, 3 N network

T = 5;      % total length of recording (seconds)
dt = 0.001; % seconds

f0 = 1/dt;  % sampling frequency (Hz)
N1 = T*f0;   % number of samples needed
df = 1/T;   % frequency resolution
fNQ = f0/2; % Nyquist frequency

taxis = dt:dt:T; % time axis
%%% simulate data ---------------------------------------------------------
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

figure;
subplot (2,3,[1 3])                              %Plot signals.
% plot(taxis,x1)
% hold on
% plot(taxis,x2)
% plot(taxis,x3)
% plotchannels([x1';x2';x3']')
ylabel('Signal')
xlabel('Time (seconds)')
legend('x1','x2','x3')
title('Sim 1','FontSize',15);

%%% Define inputs for model building and GoF -------------------------
data = [x1';x2';x3'];                      % Conglomerate data in one matrix
nlags = 40;                               % Define order of AR model
                                           % needs to be larger than true
                                           % order
                                           
data = lsfilter(data', f0, [4 50]);
figure; plotchannels(data);
data = data';
%%%% true coefficients 
true_coeffs = cell(3);
true_coeffs{1,1} = padarray(a1,nlags-L,'post'); % True coefficents up to 
true_coeffs{1,2} = padarray(a2,nlags-L,'post'); %     number of lags being estimated
true_coeffs{1,3} = zeros(nlags,1);

true_coeffs{2,1} = zeros(nlags,1);
true_coeffs{2,2} = padarray(a2,nlags-L,'post');
true_coeffs{2,3} = zeros(nlags,1);

true_coeffs{3,1} = zeros(nlags,1);
true_coeffs{3,2} = zeros(nlags,1);
true_coeffs{3,3} = padarray(a3,nlags-L,'post');

tic
[ adj_spline ] = build_ar_splines( data, nlags); % Build network using splines
splinetime  = toc;
tic
[ adj_stand] = build_ar( data, nlags );          % Build network using standard ar
standardtime = toc;



subplot 234; 
adj_true = [1 1 0; 0 1 0; 0 0 1];          % Plot true network
plotNetwork(adj_true)
title('Truth');
subplot 235;                               % Plot network using standard regressors
str = {num2str(standardtime), ' seconds' };
plotNetwork(adj_stand)
title(strcat({'Standard, '},num2str(standardtime),{' s'}))
subplot 236;                               % Plot network using spline regressors
plotNetwork(adj_spline)
title(strcat({'Spline, '},num2str(splinetime),{' s'}))


[b_est H P] = ar_gof(adj_spline, adj_stand,data, nlags,true_coeffs,taxis,2); % Run goodness of fit tests

%% Simulation 2, 3 N network
T = 5;      % total length of recording (seconds)
dt = 0.001; % seconds

f0 = 1/dt;  % sampling frequency (Hz)
N1 = T*f0;   % number of samples needed
df = 1/T;   % frequency resolution
fNQ = f0/2; % Nyquist frequency

taxis = dt:dt:T; % time axis
%%% Simulate signals --------------------------------------------------
a1 = 0.07*[hann(20)', -0.5*ones(20,1)']';   %AR coeffictients for signal 1
a2 = 0.03*[-0.5*ones(20,1)', hann(20)']';   %                  ...signal 2
a3 = -.3*ones(size(a1));                    %                  ...signal 3

L = length(a1);                             % Number of AR terms.

N = N1+L;                                 % Number of time steps.
x1 = zeros(N,1);
x2 = zeros(N,1);    
x3 = zeros(N,1);
noise = 7;                                  % noise level


                                            % Create AR model
for k=L+1:N                                 %For each time step,
    x1(k) = sum(a1.*x1(k-1:-1:k-L)) + sum(a2.*x2(k-1:-1:k-L)) +  noise*randn();
    x2(k) = sum(a2.*x2(k-1:-1:k-L)) + sum(a1.*x1(k-1:-1:k-L)) + noise*randn();
    x3(k) = sum(a3.*x3(k-1:-1:k-L)) + noise*randn();  %...x(now) depends on past.  
end
x1 = x1(L+1:end);       
x2 = x2(L+1:end);                           %Drop the first terms (0s).
x3 = x3(L+1:end);

figure(1);
subplot(2,3,[1 3])                          %Plot signals.
plot(taxis,x1)
hold on
plot(taxis,x2)
plot(taxis,x3)
ylabel('Signal')
xlabel('Time (seconds)')
legend('x1','x2','x3')

title('Scenario 2: 1 drives 1 and 2, 2 drives 1 and 2, and 3 is independent');

%%% Define inputs for model building and GoF -------------------------
data = [x1';x2';x3'];                      % Conglomerate data in one matrix
nlags = 100;                                % Define order of AR model

%%%% true coefficients

true_coeffs = cell(3);
true_coeffs{1,1} = padarray(a1,nlags-L,'post'); % True coefficents up to 
true_coeffs{1,2} = padarray(a2,nlags-L,'post'); %     number of lags being estimated
true_coeffs{1,3} = zeros(nlags,1);

true_coeffs{2,1} = padarray(a1,nlags-L,'post');
true_coeffs{2,2} = padarray(a2,nlags-L,'post');
true_coeffs{2,3} = zeros(nlags,1);

true_coeffs{3,1} = zeros(nlags,1);
true_coeffs{3,2} = zeros(nlags,1);
true_coeffs{3,3} = padarray(a3,nlags-L,'post');


tic
[ adj_spline ] = build_ar_splines( data, nlags); % Build network using splines
splinetime  = toc;
tic
[ adj_stand] = build_ar( data, nlags );          % Build network using standard AR
standardtime = toc; 

subplot 234; 
adj_true = [1 1 0; 1 1 0; 0 0 1];          % Plot true network
plotNetwork(adj_true)
title('Truth');
subplot 235;                               % Plot network using standard regressors
plotNetwork(adj_stand)
title(strcat({'Standard, '},num2str(standardtime),{' s'}))
subplot 236;                               % Plot network using spline regressors
plotNetwork(adj_spline)
title(strcat({'Spline, '},num2str(splinetime),{' s'}))

 
ar_gof(adj_spline, adj_stand,data, nlags,true_coeffs,taxis,1); % Run goodness of fit tests

%% Simulation 3, 3 N network
T = 5;      % total length of recording (seconds)
dt = 0.001; % seconds

f0 = 1/dt;  % sampling frequency (Hz)
N1 = T*f0;   % number of samples needed
df = 1/T;   % frequency resolution
fNQ = f0/2; % Nyquist frequency

taxis = dt:dt:T; % time axis

%%% Simulate signals ----------------------------------------------------
a1 = 0.07*[hann(20)', -0.5*ones(20,1)']';   % AR coefficients for signal 1
a3 = 0.03*[-0.5*ones(20,1)', hann(20)']';%-.3*ones(size(a1));       %  ... signal 3

L = length(a1);                            % Number of AR terms.

N = N1+L;                                % Number of time steps.


noise = 5 ;                                % noise level


x1 = zeros(N,1);                           % Create the AR model, 
x2 = zeros(N,1);                           
x3 = zeros(N,1);                             
                                 
for k=L+1:N 
                                                      % for each time step,
    x1(k) = sum(a1.*x3(k-1:-1:k-L)) + noise*randn();
    x2(k) = sum(a1.*x3(k-1:-1:k-L)) + noise*randn();
    x3(k) = sum(a3.*x3(k-1:-1:k-L)) + noise*randn();  %...x(now) depends on past.                                        %...x(now) depends on past.

end
x1 = x1(L+1:end);                                     % Drop the first terms (0s).
x2 = x2(L+1:end);
x3 = x3(L+1:end);                           



figure;
subplot(2,3,[1 3])                              % Plot signals.
plot(taxis,x1)
hold on
plot(taxis,x2)
plot(taxis,x3)
ylabel('Signal')
xlabel('Time (seconds)')
legend('x1','x2','x3')
title('Scenario 3: 3 drives 1 and 2');

%%% Define inputs for model building with splines -------------------------
data = [x1';x2';x3'];                      % Conglomerate data in one matrix
nlags = 40;                                % Define order of AR model

%%%% true coefficients

true_coeffs = cell(3);
true_coeffs{1,1} = zeros(nlags,1);       % True coefficents up to 
true_coeffs{1,2} = zeros(nlags,1);       %     number of lags being estimated
true_coeffs{1,3} = padarray(a1,nlags-L,'post');

true_coeffs{2,1} = zeros(nlags,1);
true_coeffs{2,2} = zeros(nlags,1);
true_coeffs{2,3} = padarray(a1,nlags-L,'post');

true_coeffs{3,1} = zeros(nlags,1);
true_coeffs{3,2} = zeros(nlags,1);
true_coeffs{3,3} = padarray(a3,nlags-L,'post');


tic
[ adj_spline ] = build_ar_splines( data, nlags); % Build network using splines
splinetime  = toc;
tic
[ adj_stand] = build_ar( data, nlags );          % Build network using standard AR
standardtime = toc;

subplot 234; 
adj_true = [0 0 1; 0 0 1; 0 0 1];          % Plot true network
plotNetwork(adj_true)
title('Truth');
          
subplot 235;                               % Plot network using standard regressors
plotNetwork(adj_stand)
title(strcat({'Standard, '},num2str(standardtime),{' s'}))

subplot 236;                               % Plot network using spline regressors
plotNetwork(adj_spline)
title(strcat({'Spline, '},num2str(splinetime),{' s'}))



ar_gof(adj_spline, adj_stand,data, nlags,true_coeffs,taxis,3); % Run goodness of fit tests

%% Simulation 4, 6 N network---------------------------------------------
T = 5;      % total length of recording (seconds)
dt = 0.001; % seconds

f0 = 1/dt;  % sampling frequency (Hz)
N1 = T*f0;   % number of samples needed
df = 1/T;   % frequency resolution
fNQ = f0/2; % Nyquist frequency

taxis = dt:dt:T; % time axis

a1 = 0.07*[0.1*hann(20)', -0.5*ones(20,1)']';
%a1 = 0.07*[hann(20)', -0.5*ones(20,1)']';   %AR coefficients for signal 1
a2 = 0.03*[-0.5*ones(20,1)', hann(20)']';    %                  ...signal 2
a3 = -0.3*ones(size(a1));                    %                  ...signal 3
a4 = 0.02*[hann(30)', -0.2*ones(1,10)]';
a5 = 0.1*[-0.8*ones(1,30), -0.8*hann(10)']';
a6 = 0.25*[-0.6*ones(20,1)', -0.2*hann(20)']'; 


L = length(a1);                             % Number of AR terms.
N = N1+L;                                 % Number of time steps.

x1 = zeros(N,1);
x2 = zeros(N,1);                           %Create the AR model,
x3 = zeros(N,1);    
x4 = zeros(N,1);
x5 = zeros(N,1);
x6 = zeros(N,1);

noise = 7;                                 % noise level

for k=L+1:N                                % For each time step,
    x1(k) = sum(a1.*x1(k-1:-1:k-L)) + sum(a2.*x2(k-1:-1:k-L)) +  noise*randn(); 
    x2(k) = sum(a1.*x2(k-1:-1:k-L)) + sum(a4.*x4(k-1:-1:k-L)) +noise*randn();
    x3(k) = sum(a3.*x3(k-1:-1:k-L)) + sum(a2.*x4(k-1:-1:k-L)) + noise*randn();  %...x(now) depends on past.

    x4(k) = sum(a4.*x4(k-1:-1:k-L)) + noise*randn();
    x5(k) = sum(a5.*x2(k-1:-1:k-L)) + sum(a4.*x5(k-1:-1:k-L)) +noise*randn();
    x6(k) = sum(a6.*x6(k-1:-1:k-L)) + sum(a4.*x3(k-1:-1:k-L)) +noise*randn();
end
x1 = x1(L+1:end);                          % Drop the first L terms (0s).   
x2 = x2(L+1:end);                             
x3 = x3(L+1:end);
x4 = x4(L+1:end);
x5 = x5(L+1:end);
x6 = x6(L+1:end);
figure;
subplot (2,3,[1 3])                              %Plot signals.
plot(taxis,x1)
hold on
plot(taxis,x2)
plot(taxis,x3)
plot(taxis,x4)
plot(taxis,x5)
plot(taxis,x6)
ylabel('Signal')
xlabel('Time')
legend('x1','x2','x3','x4','x5','x6')
title('Six Node Network','FontSize',15);

data = [x1';x2';x3';x4';x5';x6'];                      % Conglomerate data in one matrix
nlags = 40;                               % Define order of AR model
                                           % needs to be larger than true
                                           % order

%%%% true coefficients 
true_coeffs = cell(6);
true_coeffs{1,1} = padarray(a1,nlags-L,'post'); % True coefficents up to 
true_coeffs{1,2} = padarray(a2,nlags-L,'post'); %     number of lags being estimated
true_coeffs{1,3} = zeros(nlags,1);
true_coeffs{1,4} = zeros(nlags,1);
true_coeffs{1,5} = zeros(nlags,1);
true_coeffs{1,6} = zeros(nlags,1);

true_coeffs{2,1} = zeros(nlags,1);
true_coeffs{2,2} = padarray(a1,nlags-L,'post');
true_coeffs{2,3} = zeros(nlags,1);
true_coeffs{2,4} = padarray(a4,nlags-L,'post');
true_coeffs{2,5} = zeros(nlags,1);
true_coeffs{2,6} =zeros(nlags,1);

true_coeffs{3,1} = zeros(nlags,1);
true_coeffs{3,2} = zeros(nlags,1);
true_coeffs{3,3} = padarray(a3,nlags-L,'post');
true_coeffs{3,4} = padarray(a2,nlags-L,'post');
true_coeffs{3,5} = zeros(nlags,1);
true_coeffs{3,6} = zeros(nlags,1);

true_coeffs{4,1} = zeros(nlags,1);
true_coeffs{4,2} = zeros(nlags,1);
true_coeffs{4,3} = zeros(nlags,1);
true_coeffs{4,4} = padarray(a4,nlags-L,'post');
true_coeffs{4,5} = zeros(nlags,1);
true_coeffs{4,6} = zeros(nlags,1);

true_coeffs{5,1} = zeros(nlags,1);
true_coeffs{5,2} = padarray(a5,nlags-L,'post');
true_coeffs{5,3} = zeros(nlags,1);
true_coeffs{5,4} = zeros(nlags,1);
true_coeffs{5,5} = padarray(a4,nlags-L,'post');
true_coeffs{5,6} = zeros(nlags,1);

true_coeffs{6,1} = zeros(nlags,1);
true_coeffs{6,2} = zeros(nlags,1);
true_coeffs{6,3} = padarray(a4,nlags-L,'post');
true_coeffs{6,4} = zeros(nlags,1);
true_coeffs{6,5} = zeros(nlags,1);
true_coeffs{6,6} = padarray(a6,nlags-L,'post');


tic
[ adj_spline ] = build_ar_splines( data, nlags); % Build network using splines
splinetime  = toc;
tic
[ adj_stand] = build_ar( data, nlags );          % Build network using standard ar
standardtime = toc;

subplot 234;                                     % Plot true network
adj_true = [1 1 0 0 0 0; 0 1 0 1 0 0; 0 0 1 1 0 0; 0 0 0 1 0 0; 0 1 0 0 1 0; 0 0 1 0 0 1];  

plotNetwork(adj_true)
title('Truth');
subplot 235;                               % Plot network using standard regressors
str = {num2str(standardtime), ' seconds' };
plotNetwork(adj_stand)
title(strcat({'Standard, '},num2str(standardtime),{' s'}))
subplot 236;                               % Plot network using spline regressors
plotNetwork(adj_spline)
title(strcat({'Spline, '},num2str(splinetime),{' s'}))

 
b_est = ar_gof(adj_spline, adj_stand,data, nlags,true_coeffs,taxis,1); % Run goodness of fit tests

% %%
 h = get(0,'children');
 for i=1:length(h)
   saveas(h(i), ['figure' num2str(i)], 'jpg');
 end
