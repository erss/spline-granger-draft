%%% simulate network


if strcmp(model_true.noise_type,'white')
    [model_true.data] = simulate_data(model_true);
    
elseif strcmp(model_true.noise_type,'pink')
    model_true.data = make_pink_noise(0.33,model_true.T*model_true.sampling_frequency,1/model_true.sampling_frequency);
elseif strcmp(model_true.noise_type,'real')   %%%%% if real data edit this
    
    load('data_unfiltered.mat')
    data = data - repmat(mean(data,2), [1,129290]);
    N = model_true.sampling_frequency * model_true.T ;
    
    sz = 9e4:9e4+N-1;
    pre_sz = 2e4:2e4+N-1;%2e4:2e4+N-1;
    
    ntwk = model_true.ntwk; 
    data_sz= data(ntwk,sz);
    data_presz= data(ntwk,pre_sz);
    
    d = designfilt('bandstopiir','FilterOrder',2, ...
               'HalfPowerFrequency1',59,'HalfPowerFrequency2',61, ...
               'DesignMethod','butter','SampleRate',500);
           
    data_sz = filtfilt(d,data_sz')';
    data_presz = filtfilt(d,data_presz')';
    if strcmp(model_true.sztype,'sz')
        model_true.data = data_sz;
    elseif strcmp(model_true.sztype,'presz')
        model_true.data = data_presz;
    else
        fprintf('no data used')
    end
else
    fprintf('no data simulated')
end