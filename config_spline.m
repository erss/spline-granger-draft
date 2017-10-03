%%% Model type ------------------------------------------------------------
model_true.noise_type = 'white'; % 'white', 'pink', 'real'

%%% Simulation parameters -------------------------------------------------

model_true.sampling_frequency = 500;
model_true.T = 2;   % time in seconds of window
model_true.noise = 0.25;
taxis = (1/model_true.sampling_frequency):(1/model_true.sampling_frequency):model_true.T;
model_true.taxis = taxis;

if strcmp(model_true.noise_type,'white')
    %load('large_network_coef.mat');
   % load('ninenode_exp_stand.mat');
    model_true.true_coefficients =three_node_sim_1;% %%%% MODIFY COEFFICIENTS HERE!
    model_true.model_coefficients = model_true.true_coefficients;   
elseif strcmp(model_true.noise_type,'real')
    model_true.sztype = 'sz'; % presz
    nwk =  1:(32+6);
    badchannels = [1,9,21,32, 8,31]; % % badchannels = [1,9,21,32,83, 8,31];
    nwk(badchannels) = [];
    
    model_true.ntwk = nwk;
 %  model_true.ntwk = [2 7 18 20 21 46 50 55 60 65 70 90 80 72 ];     
end
%%% Define model inputs for spline Granger & standard Granger -------------

model_true.s = 0.5;                     % tension parameter for spline
model_true.estimated_model_order = 20;  % model_order used to estimate

%number_of_knots      = floor(model_true.estimated_model_order/3);
model_true.cntrl_pts =[0:5:model_true.estimated_model_order]; %make_knots(model_true.estimated_model_order,number_of_knots);

%%% Define network testing parameters -------------------------------------

model_true.q = 0.05;            % FDR max number acceptable proportion of false discoveries
model_true.nsurrogates = 100;   % number of surrogates used for bootstrapping
model_true.nrealizations = 3; % number of realizations used for spectral testing


