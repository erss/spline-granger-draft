%%% Generate data -------------------------------------------------------
T = 5;      % total length of recording (seconds)
dt = 0.001; % seconds

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
%%% Build network using splines -----------------------------------------

[ adj_spline ] = build_ar_splines( data, nlags); 
b = estimate_coef(data,adj_spline, nlags); % estimated coefficients  from adj_spline

%%% Generate surrogate data given adj_spline ----------------------------

which_electrode = 1; % which electrode in network to generate data for 

nsurrogates = 1000; % number of surrogates

b_trial = generate_surrogates(data,adj_spline,nlags,which_electrode,nsurrogates);

figure;
for row = 1:nsurrogates
 plot(b_trial(row,:))  % plot surrogate data
 hold on;
end

plot(b(which_electrode,:),'--k','LineWidth',1)

%%% Bootstrap for each coefficient estimate on the surrogate data
bounds = zeros(2,size(b_trial,2));
nshuffles = 10000; % number of times you resample from data
nsamples = nsurrogates*0.75; % number of samples you take
for k = 1:size(b_trial,2)
    
    bounds(:,k) = myBootstrap(b_trial(:,k),nshuffles, nsamples);
end


%%% plot estimated coefficients and confidence bounds
figure;
plot(b(which_electrode,:),'r');
hold on
plot((bounds(1,:)),'--b');
plot((bounds(2,:)),'--r');
plot([a1' a2' zeros(1,40)])


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%% Simulation 1, 3 N network
% 
% 
% T = 5;      % total length of recording (seconds)
% dt = 0.001; % seconds
% 
% f0 = 1/dt;  % sampling frequency (Hz)
% N1 = T*f0;   % number of samples needed
% df = 1/T;   % frequency resolution
% fNQ = f0/2; % Nyquist frequency
% 
% taxis = dt:dt:T; % time axis
% 
% nlags = 40;                               % Define order of AR model
%                                            % needs to be larger than true
%                                            % order
% ntrials = 100;
% nelectrodes = 6;
% b_trial = zeros(ntrials,nlags*nelectrodes);
% tester = zeros(ntrials,nelectrodes);
% 
% 
% a1 = 0.07*[hann(20)', -0.5*ones(20,1)']';   %AR coefficients for signal 1
% a2 = 0.03*[-0.5*ones(20,1)', hann(20)']';    %                  ...signal 2
% a3 = -0.3*ones(size(a1));                    %                  ...signal 3
% a4 = 0.02*[hann(30)', -0.2*ones(1,10)]';
% a5 = 0.1*[-0.8*ones(1,30), -0.8*hann(10)']';
% a6 = 0.25*[-0.6*ones(20,1)', -0.2*hann(20)']'; 
% 
% 
% L = length(a1);                             % Number of AR terms.
% N = N1+L;                                 % Number of time steps.
% 
% x1 = zeros(N,1);
% x2 = zeros(N,1);                           %Create the AR model,
% x3 = zeros(N,1);    
% x4 = zeros(N,1);
% x5 = zeros(N,1);
% x6 = zeros(N,1);
% 
% noise = 10;                                 % noise level
% 
% 
% for trial = 1:ntrials
%     for k=L+1:N                                % For each time step,
%         x1(k) = sum(a1.*x1(k-1:-1:k-L)) + sum(a2.*x2(k-1:-1:k-L)) +  noise*randn(); 
%         x2(k) = sum(a1.*x2(k-1:-1:k-L)) + sum(a4.*x4(k-1:-1:k-L)) +noise*randn();
%         x3(k) = sum(a3.*x3(k-1:-1:k-L)) + sum(a2.*x4(k-1:-1:k-L)) + noise*randn();  %...x(now) depends on past.
% 
%         x4(k) = sum(a4.*x4(k-1:-1:k-L)) + noise*randn();
%         x5(k) = sum(a5.*x2(k-1:-1:k-L)) + sum(a4.*x5(k-1:-1:k-L)) +noise*randn();
%         x6(k) = sum(a6.*x6(k-1:-1:k-L)) + sum(a4.*x3(k-1:-1:k-L)) +noise*randn();
%     end
%     
%     x1 = x1(L+1:end);                          % Drop the first L terms (0s).   
%     x2 = x2(L+1:end);                             
%     x3 = x3(L+1:end);
%     x4 = x4(L+1:end);
%     x5 = x5(L+1:end);
%     x6 = x6(L+1:end);
%     
%     data = [x1';x2';x3';x4';x5';x6'];         % Put data in one matrix
%     
%    [adj_mat] = build_ar_splines( data, nlags); % Use F-test to determine connectivity
%     
%    %%% Determine coefficients assuming adj_mat connectivity
%    
%     
%     nobservations = length(data(1,nlags+1:end)); % number of observations
% 
% 
%     %%% Define control points and build predictors
% 
%     c_pt_times = [0:10:nlags];  % Define Control Point Locations
% 
% 
%     s = 0.5;                                    % Define Tension Parameter
% 
%     % Construct spline regressors for case nelectrodes = 1.
%     c_pt_times_all = [c_pt_times(1)-2 c_pt_times c_pt_times(end)+2];
%     Z = zeros(nlags,length(c_pt_times_all));
%     num_c_pts = length(c_pt_times_all);  %number of control points in total
%     for i=1:nlags
%         nearest_c_pt_index = max(find(c_pt_times_all<i));
%         nearest_c_pt_time = c_pt_times_all(nearest_c_pt_index);
%         next_c_pt_time = c_pt_times_all(nearest_c_pt_index+1);
%         u = (i-nearest_c_pt_time)/(next_c_pt_time-nearest_c_pt_time);
%         p=[u^3 u^2 u 1]*[-s 2-s s-2 s;2*s s-3 3-2*s -s;-s 0 s 0;0 1 0 0];
%         Z(i,:) = [zeros(1,nearest_c_pt_index-2) p zeros(1,num_c_pts-4-(nearest_c_pt_index-2))];
%     end
% 
% 
% 
% %%% Build model for only electrode ii in network
%    ii = 2;
%    
%    
%     % Build history regressors
%     preds = logical(adj_mat(ii,:)); %% Use results of F-test to input in network
%     data_copy = data(preds,:); %% Remove electrodes not connected in spline network
%     
%     
%     Z1 = kron(eye(size(data_copy,1)),Z);     % Build block diagonal spline regressors
%   
%      X = [];                                 % Build history matrix
%     for k = 1:size(data_copy,1)
%         X_temp = []; 
%         sgnl = data_copy(k,:)';
%         for i=1:nlags                                   %For each lag,
%             X_temp = [X_temp, circshift(sgnl,i)];   %... shift x and store it.
%         end
%         X_temp = X_temp(nlags+1:end,:);  
%         X = [X X_temp];
%     end
%     
%     
%     %%% Build Model
%      
%     % Generate observations for given y
%     y = data(ii,nlags+1:end);   
%     y = y';
% 
%     % Fit full model and calculate RSS
%     Xfull = X * Z1;      % regressors for y_hat = X*Z1*alpha
%     [alpha,~,stats] = glmfit(Xfull,y);  % estimate values at control points, alpha
%     bhat = Z1*alpha(2:end);             % calculate beta values, for every point in space
%                                        % only for electrodes in network
% 
% 
% 
%     j =1;
%     for k = 1:nelectrodes
% 
%        if preds(k)
%        b_trial(trial,((k-1)*nlags + 1): k*nlags) = bhat(((j-1)*nlags + 1): j*nlags);         
% 
%        j = j+1;
%        end
% 
%     end
%     tester(trial,:) = preds;
%     
%     
%     
% end

