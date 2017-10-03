%%%% Real Data!
clear all;
%%% Model type ------------------------------------------------------------
model_true.noise_type = 'real'; % 'white', 'pink', 'real'

 model_true.sztype = 'presz'; % presz
 model_true.estimated_model_order = 20;  % model_order used to estimate


%%% Simulation parameters -------------------------------------------------

model_true.sampling_frequency = 500;
model_true.T = 2;   % time in seconds of window
model_true.noise = 0.25;
taxis = (1/model_true.sampling_frequency):(1/model_true.sampling_frequency):model_true.T;
model_true.taxis = taxis;

%%% Define model inputs for spline Granger & standard Granger -------------

model_true.s = 0.5;                     % tension parameter for spline
model_true.estimated_model_order = 20;  % model_order used to estimate

number_of_knots      = floor(model_true.estimated_model_order/3);
model_true.cntrl_pts = [0 5:5:model_true.estimated_model_order];%make_knots(model_true.estimated_model_order,number_of_knots);

%%% Define network testing parameters -------------------------------------

model_true.q = 0.05;            % FDR max number acceptable proportion of false discoveries
model_true.nsurrogates = 1000;   % number of surrogates used for bootstrapping
model_true.nrealizations = 20; % number of realizations used for spectral testing

 badchannels = [1,8,9,13,20,21,22,23,24,25,26,31,32,34,38,41,48,49,50,68,69,71,77,78,79,82,83,87,88 ,89];
 ntwork =     1:94;
 ntwork(badchannels)=[];

for k = 53 %1:64;%[20 49 53 54]%1:64
 model_true.ntwk =ntwork(k);%76;
 simulate_network;
 infer_network;

figure
subplot 221

plot(model_true.taxis,model_true.data,'k','LineWidth',2)
     hold on;
     plot([0,.2], [min(model_true.data(1,:)), min(model_true.data(1,:))], 'k', 'LineWidth', 2);
set(gca,'YTickLabel',[]);
set(gca,'XTickLabel',[]);
set(gca,'XTick',[]);
set(gca,'YTick',[]);

title('Single Electrode Recording','FontSize',20);
box off


ylabel('Signal','FontSize',18);
xlabel('Time (s)','FontSize',18);

%title(num2str(i),'FontSize',20);

subplot 222

list=lines(size(model_true.data,1));
mySpec(model_true.data(1,:),model_true.sampling_frequency,'yesplot','notapers','k');
box off
ylim([-60 90])
% a = get(gca,'YTickLabel');
% set(gca,'YTickLabel',a,'fontsize',16)
title('Spectrogram','FontSize',20);

% 
 subplot(2,2,[3 4])
 gof_bootstrap(model_true,model_spline,model_standard);
box off
% a = get(gca,'YTickLabel');
% set(gca,'YTickLabel',a,'fontsize',16)

title('Estimated Coefficient Fits','FontSize',20);
plot([0 0.04],[0 0],'color',[.57 .57 .57],'LineWidth',1.7)
end
%title(strcat({'Spline, '},num2str(model_spline.computation_time),{' s'},' Overlap, ',num2str(model_spline.accuracy)))
% bhat = model_standard.model_coefficients;
% % save(num2str(ntwork(i)),'bhat')
   [notwhitestand, dwstand]=dwstat( model_standard);
   [notwhitespline,dwspline]=dwstat( model_spline);


% 
h = get(0,'children');
for i=1:length(h)
saveas(h(i), ['fig5_single'  num2str(i)], 'fig');     
end
close all;

save('fig5_singlenode')