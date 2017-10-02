%%%%
clear all;
ntrials=2;

dw_fails_spline =[];
dw_fails_standard =[];
gr_fails_spline =[];
gr_fails_standard =[];

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
%   load('/Users/erss/Documents/MATLAB/ar_model/Simulation_coefficients/nine_node/b_standard_order35_rdi.mat');
load('bhat_stand30_tues.mat');
model_true.true_coefficients =bhat;
%nine_node_order20_rdi; %%%% MODIFY COEFFICIENTS HERE!
    model_true.model_coefficients = model_true.true_coefficients;   
end
%%% Define model inputs for spline Granger & standard Granger -------------

model_true.s = 0.5;                     % tension parameter for spline
model_true.estimated_model_order = 30;  % model_order used to estimate
model_true.cntrl_pts = [0:5:model_true.estimated_model_order];%make_knots(model_true.estimated_model_order,number_of_knots);

%%% Define network testing parameters -------------------------------------

model_true.q = 0.05;            % FDR max number acceptable proportion of false discoveries
model_true.nsurrogates = 1000;   % number of surrogates used for bootstrapping
model_true.nrealizations = 20; % number of realizations used for spectral testing


nelectrodes = size(model_true.model_coefficients,1);

trials_spline = zeros(nelectrodes,nelectrodes,ntrials); 
trials_stand = zeros(nelectrodes,nelectrodes,ntrials);

splinetimes = zeros(1,ntrials);
standardtimes = zeros(1,ntrials);

for k = 1:ntrials
    
    simulate_network;
    
%     tic
%     [ adj_spline] = build_ar_splines( model_true);
%     splinetime(k)  = toc;
%     
%     tic
%     [ adj_mat] = build_ar( model_true);
%     standardtime(k)  = toc;
%     
%     trials_spline(:,:,k) = adj_spline;
%     trials_stand(:,:,k) = adj_mat;

tic
[ adj_spline] = build_ar_splines( model_true);
splinetime  = toc;
[ bhat, yhat] = estimate_coefficient_fits( model_true, adj_spline);


model_spline = model_true;
model_spline.model_coefficients = bhat;
model_spline.computation_time = splinetime;
model_spline.signal_estimate = yhat;
model_spline.network = adj_spline;


     trials_spline(:,:,k) = model_spline.network;
     
     if k == 1
        tic
        [ adj_mat] = build_ar( model_true);
        standardtime(k)  = toc;
        trials_stand(:,:,1) = adj_mat;
     end
     splinetimes(k)= model_spline.computation_time;
     %standardtimes(k)=model_standard.computation_time;

     
         [nw] = dwstat(model_spline);
    dw_fails_spline = [dw_fails_spline nw];
     
%          [nw] = dwstat(model_standard);
%     dw_fails_standard = [dw_fails_standard nw];

%     clear m2fit
%     clear m3fit
%     [ m2fit,m3fit] = grstat1( model_true,model_standard,model_spline );
%     gr_fails_standard = [ gr_fails_standard m2fit.fails'];
%     gr_fails_spline = [ gr_fails_spline m3fit.fails'];
    
%       [m2fit, m3fit] = grstat1(model_true,model_spline,model_standard);
%  ts_spline(i) = m2fit.stat;
%  ts_stand(i)  = m3fit.stat;
 
    fprintf([num2str(k) ' \n']);
end

avg_network = mean(trials_spline,3);
std_network = std(trials_spline,0,3);
b = model_true.model_coefficients;
%%
%%% Figure 1
figure;
%%% plot coefficient strengths
subplot 221
[bb, ii] = max(abs(b),[],3);
plotNetwork(bb)
title('Coefficient Strength','FontSize',20)
colorbar
caxis([0 1])

%%% plot one instance of spline inference
subplot 222
plotNetwork(trials_spline(:,:,1))
title('Spline Network','FontSize',20)
colorbar
caxis([0 1])

%%% plot average network
subplot 223
plotNetwork(avg_network);
title('Averaged Network','FontSize',20);
colorbar
caxis([0 1])

%%% plot std network
subplot 224
plotNetwork(std_network)
title('Standard Deviation','FontSize',20);
colorbar
caxis([0 1])
%%
adj_true = model_true.true_coefficients;
adj_true(adj_true~=0)=1;
adj_true=sum(adj_true,3);
adj_true(adj_true~=0)=1;

adj_thresh = bb;
adj_thresh(adj_thresh >= 0.1)= 1;
adj_thresh(adj_thresh < 0.1)= 0;

for i = 1:ntrials
    accspline(i) = network_accuracy(adj_true,trials_spline(:,:,i));
 %  accstand(i) = network_accuracy(adj_true,trials_stand(:,:,i));
    
    threshspline(i) = network_accuracy(adj_thresh,trials_spline(:,:,i));
    threshstand(i) = network_accuracy(adj_thresh,trials_stand(:,:,i));
    
end
% 
% %%% Figure 2
% Labels = {'Standard', 'Spline'};
% figure;
% %%% computational time bar plot
% subplot 121
% barplot(Labels, standardtimes(2:end),splinetimes(2:end))
% set(gca,'xlim',[0.5 2.5])
% ylabel('Computation time (s)','FontSize',18)
% 
% %%% accuracy bar plot
% subplot 122
% barplot(Labels,accstand,accspline)
% ylabel('Accuracy','FontSize',18)
% set(gca,'xlim',[0.5 2.5])
% ylim([0 1])

%%% threshold accuracy bar plot
% subplot 133
% barplot(Labels,threshstand,threshspline)
% ylabel('Accuracy (Threshold)','FontSize',18)
% set(gca,'xlim',[0.5 2.5])
% ylim([0 1])
%% Figure 3
figure;

%%% plot signals
subplot 211
plotchannels(model_true.taxis,model_true.data')

ylabel('Signal','FontSize',18);
xlabel('Time (s)','FontSize',18);
title('Nine Node Network Simulation','FontSize',20);

%%% plot spectrogram
subplot 212
list=lines(9);
for i =1:9
    mySpec(model_true.data(i,:),model_true.sampling_frequency,'yesplot','tapers',list(i,:));
    hold on;
    
end
%ylim([-100,0]);
title('Spectrogram','FontSize',20);
ylabel('Power (dB)','FontSize',18)
xlabel('Frequency (Hz)','FontSize',18)
box off
% a = get(gca,'YTickLabel');
% set(gca,'YTickLabel',a,'fontsize',16)


%% SAVE ALL ----------------------------
% save('fig3_standard')
% 
% h = get(0,'children');
% for i=1:length(h)
%     
%     saveas(h(i), ['fig3_avgntwk_standardrdi'  num2str(i)], 'fig');   
% end
% close all;


