%%% Model type ------------------------------------------------------------
model_true.noise_type = 'real'; % 'white', 'pink', 'real'

%%% Simulation parameters -------------------------------------------------

model_true.sampling_frequency = 500;
model_true.T = 2;   % time in seconds of window
model_true.noise = 0.25;
taxis = (1/model_true.sampling_frequency):(1/model_true.sampling_frequency):model_true.T;
model_true.taxis = taxis;

if strcmp(model_true.noise_type,'white')
    %load('large_network_coef.mat');
    model_true.true_coefficients =single_node_order20;% %%%% MODIFY COEFFICIENTS HERE!
    model_true.model_coefficients = model_true.true_coefficients;   
elseif strcmp(model_true.noise_type,'real')
    model_true.sztype = 'sz'; % presz
    nwk =  1:(32+6);
    badchannels = [1,9,21,32, 8,31]; % % badchannels = [1,9,21,32,83, 8,31];
    nwk(badchannels) = [];
    
    model_true.ntwk = nwk;
    model_true.ntwk = [2 7 18 22 42 46 90 80 77 ];     
end
%%% Define model inputs for spline Granger & standard Granger -------------

model_true.s = 0.5;                     % tension parameter for spline
model_true.estimated_model_order = 35;  % model_order used to estimate

number_of_knots      = floor(model_true.estimated_model_order/3);
model_true.cntrl_pts = make_knots(model_true.estimated_model_order,number_of_knots);

%%% Define network testing parameters -------------------------------------

model_true.q = 0.05;            % FDR max number acceptable proportion of false discoveries
model_true.nsurrogates = 15;   % number of surrogates used for bootstrapping
model_true.nrealizations = 20; % number of realizations used for spectral testing


