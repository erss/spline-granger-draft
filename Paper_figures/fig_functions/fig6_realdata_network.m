%%%% Real Data!
clear all;
%%% Model type ------------------------------------------------------------
model_true.noise_type = 'real'; % 'white', 'pink', 'real'

%%% Simulation parameters -------------------------------------------------

model_true.sampling_frequency = 500;
model_true.T = 2;   % time in seconds of window
model_true.noise = 0.25;
taxis = (1/model_true.sampling_frequency):(1/model_true.sampling_frequency):model_true.T;
model_true.taxis = taxis;

if strcmp(model_true.noise_type,'real')
    model_true.sztype = 'sz'; % presz
     badchannels = [1,8,9,13,20,21,22,23,24,25,26,31,32,34,38,41,48,49,50,68,69,71,77,78,79,82,83,87,88,89];
 ntwk =     1:94;
 ii = [ 57,63,56,50,39,45,26,54,17];
 ntwk(badchannels)=[];
 ntwk=ntwk(ii);
     model_true.ntwk = ntwk;%[2 7 18 22 42 46 90 80 77];     % badchannels = [1,9,21,32,83, 8,31];

 %  model_true.ntwk = [2 7 18 20 21 46 50 55 60 65 70 90 80 72 ];     
end
%%% Define model inputs for spline Granger & standard Granger -------------

model_true.s = 0.5;                     % tension parameter for spline
model_true.estimated_model_order = 20;  % model_order used to estimate

number_of_knots      = floor(model_true.estimated_model_order/3);
model_true.cntrl_pts = make_knots(model_true.estimated_model_order,number_of_knots);

%%% Define network testing parameters -------------------------------------

model_true.q = 0.05;            % FDR max number acceptable proportion of false discoveries
model_true.nsurrogates = 1000;   % number of surrogates used for bootstrapping
model_true.nrealizations = 10; % number of realizations used for spectral testing



 
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

figure;
subplot 121
plotNetwork(model_standard.network);
title('Standard Granger','FontSize',20)
disp([num2str(model_standard.computation_time) '']);
subplot 122
plotNetwork(model_spline.network);
title('Spline Granger','FontSize',20)
disp([num2str(model_spline.computation_time) '']);
disp(num2str(model_spline.accuracy));

%title(strcat({'Spline, '},num2str(model_spline.computation_time),{' s'},' Overlap, ',num2str(model_spline.accuracy)))
nwspline= dwstat( model_spline);
nwstandard =dwstat( model_standard);

%%% SAVE ALL ----------------------------
save('fig6')


h = get(0,'children');
for i=1:length(h)
saveas(h(i), ['fig6_ntwk'  num2str(i)], 'fig');     
end
close all;