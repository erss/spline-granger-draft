function data = simulate_data(model)
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

    for k = nlags:length(data)-1
        data(:,k+1) = myPrediction(data(:,1:k),b);
        data(:,k+1) = data(:,k+1) + noise.*randn(size(data,1),1);
    end

    data = data(:,nlags+3001:end);
end
