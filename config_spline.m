%%% Model type
model_true.noise_type = 'white'; % 'white', 'real'

%%% Simulation parameters

model_true.sampling_frequency = 500;
model_true.T = 2;   % time in seconds of window
model_true.noise = 0.25; 

if strcmp(model_true.noise_type,'white')
    model_true.true_coefficients = three_node_sim_1; %%%% MODIFY COEFFICIENTS HERE!
    model_true.model_coefficients = model_true.true_coefficients;
end
%%% Define model inputs for spline Granger & standard Granger

model_true.s = 0.5;                     % tension parameter for spline
model_true.estimated_model_order = 50;  % model_order used to estimate

number_of_knots       = floor(model_true.estimated_model_order/3);
model_true.cntrl_pts = make_knots(model_true.estimated_model_order,number_of_knots);

%%% Define network testing parameters

model_true.q = 0.05; % FDR max number acceptable proportion of false discoveries

%%% Define inputs for model testing
model_true.nsurrogates = 100; % number of surrogates used for bootstrapping
model_true.nrealizations = 100; % number of realizations used for spectral testing


