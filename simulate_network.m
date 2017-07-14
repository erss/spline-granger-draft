%%% simulate network


if strcmp(model_true.noise_type,'white')
    [model_true.data] = simulate_data(model_true);
   
elseif strcmp(model_true.noise_type,'pink')
   model_true.data = make_pink_noise(0.33,model_true.T*model_true.sampling_frequency,1/model_true.sampling_frequency);
elseif strcmp(model_true.noise_type,'real')   %%%%% if real data edit this
    load('data_unfiltered.mat')
%load('data.mat')
data = data - repmat(mean(data,2), [1,129290]);
data_type = 'real';

N = model_true.sampling_frequency * model_true.T ;
sz = 9e4:9e4+N-1;
pre_sz = 2e4:2e4+N-1;
ntwk=34;
% badchannels = [1,9,21,32,83, 8,31];
ntwk = [2 6 18 22 42 46 90 82 77 ];
data_sz= data(ntwk,sz);
data_presz= data(ntwk,pre_sz);

model_true.data = data_sz;
else
    fprintf('no data simulated')
end