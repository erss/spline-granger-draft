%%%% Real Data!
 
clear all;

i=32;
%%%% Ten seconds of data, compute both spline and standard ----------------
 %%% Model type ------------------------------------------------------------
model_true.noise_type = 'real'; % 'white', 'pink', 'real'
model_true.sztype = 'presz'; % presz

%%% Simulation parameters -------------------------------------------------

model_true.sampling_frequency = 500;
model_true.T = 2;   % time in seconds of window
model_true.noise = 0.25;

%%% Define model inputs for spline Granger & standard Granger -------------

model_true.s = 0.5;                     % tension parameter for spline
model_true.estimated_model_order = 50;  % model_order used to estimate


model_true.cntrl_pts = [0:5:model_true.estimated_model_order];
%%% Define network testing parameters -------------------------------------

model_true.q = 0.05;            % FDR max number acceptable proportion of false discoveries
model_true.nsurrogates = 1000;   % number of surrogates used for bootstrapping
model_true.nrealizations = 20; % number of realizations used for spectral testing


 badchannels = [1,8,9,13,20,21,22,23,24,25,26,31,32,34,38,41,48,49,50,68,69,71,77,78,79,82,83,87,88,89];
 ntwk =     1:94;
 ntwk(badchannels)=[];
 model_true.ntwk = ntwk(1:i);
 
 simulate_network;
 %infer_network;
 
 %%% Fit spline to data ---------------------------------------------------
tic
[ adj_spline] = build_ar_splines( model_true);
splinetime  = toc;

model_spline = model_true;
model_spline.computation_time = splinetime;
model_spline.network = adj_spline;
figure;
subplot 122
plotNetwork(model_spline.network);
title('Spline Granger - 2 s','FontSize',20)

%%% INFER SPLINE NETWORK FOR 2 SECONDS OF DATA ----------------------------

clear model_true;

model_true.noise_type = 'real'; % 'white', 'pink', 'real'
model_true.sztype = 'presz'; % presz

%%% Simulation parameters -------------------------------------------------

model_true.sampling_frequency = 500;
model_true.T = 10;   % time in seconds of window
model_true.noise = 0.25;

%%% Define model inputs for spline Granger & standard Granger -------------

model_true.s = 0.5;                     % tension parameter for spline
model_true.estimated_model_order = 50;  % model_order used to estimate


model_true.cntrl_pts = [0:5:model_true.estimated_model_order];
%%% Define network testing parameters -------------------------------------

model_true.q = 0.05;            % FDR max number acceptable proportion of false discoveries
model_true.nsurrogates = 1000;   % number of surrogates used for bootstrapping
model_true.nrealizations = 20; % number of realizations used for spectral testing


 badchannels = [1,8,9,13,20,21,22,23,24,25,26,31,32,34,38,41,48,49,50,68,69,71,77,78,79,82,83,87,88,89];
 ntwk =     1:94;
 ntwk(badchannels)=[];
 model_true.ntwk = ntwk(1:i);
 
 simulate_network;
 %infer_network;
 
 %%% Fit spline to data ---------------------------------------------------
tic
[ adj_spline2] = build_ar_splines( model_true);
splinetime  = toc;


subplot 121
plotNetwork(adj_spline2);
title('Spline Granger - 10 s','FontSize',20)

%%
figure;
subplot 121
plotNetwork(model_spline_two.network);
title('Spline Granger - 2 s','FontSize',20)
%subplot 132
plotNetwork(model_standard.network);
title('Standard Granger - 10 s','FontSize',20)

subplot 122
plotNetwork(model_spline.network);
title('Spline Granger - 10 s','FontSize',20)

%% calc of overlap

i=32;
indices = 1:i+1:i^2;
splinenet = adj_spline2;
standnet = adj_stand2;
splinenet(indices)=nan;
standnet(indices)=nan;
splinenet=splinenet(:);
standnet=standnet(:);

splinenet(isnan(splinenet))=[];
standnet(isnan(standnet))=[];
diff = splinenet-standnet;
frac =  length(find(diff==0))/length(diff)
%%
%nw_spline_two = dwstat( model_spline_two);


% %%% ---- SAVE IT ALL
h = get(0,'children');
for i=1:length(h)
saveas(h(i), ['fig5'  num2str(i)], 'fig');     
end
close all;

save('fig5')
beep