%%%% Real Data!
 
clear all;
%%%% Ten seconds of data, compute both spline and standard ----------------
 %%% Model type ------------------------------------------------------------
model_true.noise_type = 'real'; % 'white', 'pink', 'real'
model_true.sztype = 'presz'; % presz

%%% Simulation parameters -------------------------------------------------

model_true.sampling_frequency = 500;
model_true.T = 10;   % time in seconds of window
model_true.noise = 0.25;
taxis = (1/model_true.sampling_frequency):(1/model_true.sampling_frequency):model_true.T;
model_true.taxis = taxis;


%%% Define model inputs for spline Granger & standard Granger -------------

model_true.s = 0.5;                     % tension parameter for spline
model_true.estimated_model_order = 20;  % model_order used to estimate


model_true.cntrl_pts = [0:5:model_true.estimated_model_order];
%%% Define network testing parameters -------------------------------------

model_true.q = 0.05;            % FDR max number acceptable proportion of false discoveries
model_true.nsurrogates = 1000;   % number of surrogates used for bootstrapping
model_true.nrealizations = 20; % number of realizations used for spectral testing


 badchannels = [1,8,9,13,20,21,22,23,24,25,26,31,32,34,38,41,48,49,50,68,69,71,77,78,79,82,83,87,88,89];
 ntwk =     1:94;
 ntwk(badchannels)=[];
 model_true.ntwk = ntwk;
 
 simulate_network;
 infer_network;
 
 figure;
subplot 211
plotchannels(model_true.taxis,model_true.data')
set(gca,'YTickLabel',[]);
set(gca,'XTickLabel',[]);
set(gca,'XTick',[]);
set(gca,'YTick',[]);
box off
ylabel('Signal','FontSize',18);
xlabel('Time (s)','FontSize',18);
title('Nine Electrode Recordings','FontSize',20);


subplot 212

list=lines(size(model_true.data,1));
for i =1:size(model_true.data,1)
    mySpec(model_true.data(i,:),model_true.sampling_frequency,'yesplot','tapers',list(i,:));
    hold on;
  
end
box off
ylim([-60 120])
a = get(gca,'YTickLabel');
set(gca,'YTickLabel',a,'fontsize',16)
title('Spectrogram','FontSize',20);

% 
% 

%title(strcat({'Spline, '},num2str(model_spline.computation_time),{' s'},' Overlap, ',num2str(model_spline.accuracy)))
nw_spline_ten = dwstat( model_spline);
nw_standard = dwstat( model_standard);

%%% INFER SPLINE NETWORK FOR 2 SECONDS OF DATA ----------------------------


model_true.T = 2;   % time in seconds of window
taxis = (1/model_true.sampling_frequency):(1/model_true.sampling_frequency):model_true.T;
model_true.taxis = taxis;

simulate_network;

tic
[ adj_spline2] = build_ar_splines( model_true);
splinetime2  = toc;
[ bhat2,yhat2] = estimate_coefficient_fits( model_true, adj_spline2);

model_spline_two = model_spline;
model_spline_two.model_coefficients = bhat2;
model_spline_two.computation_time = splinetime2;
model_spline_two.signal_estimate = yhat2;
model_spline_two.network = adj_spline2;

%%
figure;
subplot 131
plotNetwork(model_spline_two.network);
title('Spline Granger - 2 s','FontSize',20)
subplot 132
plotNetwork(model_standard.network);
title('Standard Granger - 10 s','FontSize',20)
disp([num2str(model_standard.computation_time) '']);
subplot 133
plotNetwork(model_spline.network);
title('Spline Granger - 10 s','FontSize',20)
disp([num2str(model_spline.computation_time) '']);
disp(num2str(model_spline.accuracy));
%%
nw_spline_two = dwstat( model_spline_two);


%%% ---- SAVE IT ALL
h = get(0,'children');
for i=1:length(h)
saveas(h(i), ['fig7_ntwk'  num2str(i)], 'fig');     
end
close all;

save('fig7')