%%% simulate network


if strcmp(model_true.noise_type,'white')
    model_true.data = simulate_data(model_true);
    
elseif strcmp(model_true.noise_type,'pink')
    %%%%% include code to simulate pink noise data
% elseif strcmp(model_true.data_type,'real')   %%%%% if real data edit this
else
    fprintf('no data simulated')
end