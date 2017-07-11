%%% simulate network


if strcmp(model_true.noise_type,'white')
    model_true.data = simulate_data(model_true);
   
elseif strcmp(model_true.noise_type,'pink')
   model_true.data = make_pink_noise(0.33,model_true.T*model_true.sampling_frequency,1/model_true.sampling_frequency);
% elseif strcmp(model_true.noise_type,'real')   %%%%% if real data edit this
else
    fprintf('no data simulated')
end