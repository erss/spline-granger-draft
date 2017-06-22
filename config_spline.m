%%% Model type
model_true.noise_type = 'white';
model_true.data_type = 'simulation';

%%% Simulation parameters

model_true.sampling_frequency = 500;
model_true.T = 2; 
model_true.noise = 0.25;
model_true.true_coefficients = three_node_sim_1; %nine_node_order50_rdi;
model_true.model_coefficients = model_true.true_coefficients; 
%%% Define model inputs for spline Granger & standard Granger

model_true.s = 0.5;                     % tension parameter for spline
model_true.estimated_model_order = 50;  % model_order used to estimate

number_of_knots       = floor(model_true.estimated_model_order/3);
model_true.cntrl_pts = make_knots(model_true.estimated_model_order,number_of_knots);

%%% Define network testing parameters

model_true.q = 0.05; % FDR max number acceptable proportion of false discoveries

%%% Define inputs for model testing

model_true.nsurrogates = 100;
model_true.nrealizations = 1000;



%%
%%%%%%%%%%%%%%%%%%%%%%%%% REAL DATA PARAM
%%%%%%%%%%%%%%%%%%%%%%%%% white v pink noise
