function [X] = simulate_data(model)
% SIMULATE_DATA simulates data using AR model.
seed=0;
%rng_seed(seed);
    b = model.model_coefficients;
    nelectrodes = size(b,1); % number of electrodes
    nlags = size(b,3);  % true model order
    T = model.T;      % total length of recording (seconds)
    f0 = model.sampling_frequency;
    dt = 1/f0; % seconds
    df = 1/T;   % frequency resolution
    fNQ = f0/2; % Nyquist frequency

    N = T*f0 + nlags + 3000;
    noise = model.noise;
    data = zeros(nelectrodes,N);
    X = noise.*randn([nelectrodes N]);
    
%     for k = nlags:length(data)-1
%         data(:,k+1) = myPrediction(data(:,1:k),b);
%         noisely(k) = noise.*randn(size(data,1),1);
%         data(:,k+1) = data(:,k+1) + noisely(k);
%     end
%X = [0 noisely];
for t = nlags+1:N
        for k = 1:nlags
            X(:,t) = X(:,t) + b(:,:,k)*X(:,t-k);
        end
end
    X = X(:,nlags+3001:end);
%    data = data(:,nlags+3001:end);
%     
%     mySpec(X,500,'yesplot','tapers'); 
%     figure; mySpec(data,500,'yesplot','tapers','r');
end
