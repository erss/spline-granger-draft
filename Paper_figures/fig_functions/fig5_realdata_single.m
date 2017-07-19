%%%% Real Data!

config_spline;
model_true.noise_type = 'real'; % 'white', 'pink', 'real'
 model_true.sztype = 'sz'; % presz
 model_true.ntwk = 22;%[2 7 18 22 42 46 90 80 77 ]; % i;%[77];     % badchannels = [1,9,21,32,83, 8,31];
%  [2 7 18 22 42 46 90 80 77 ]; 
 simulate_network;
 infer_network;

 
 figure;
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
a = get(gca,'YTickLabel');
set(gca,'YTickLabel',a,'fontsize',16)
title('Spectrogram','FontSize',20);

% 
 subplot(2,2,[3 4])
 gof_bootstrap(model_true,model_spline,model_standard);
box off
a = get(gca,'YTickLabel');
set(gca,'YTickLabel',a,'fontsize',16)
title('Estimated Coefficient Fits','FontSize',20);


%title(strcat({'Spline, '},num2str(model_spline.computation_time),{' s'},' Overlap, ',num2str(model_spline.accuracy)))
h = get(0,'children');
for i=1:length(h)
saveas(h(i), ['fig5_single'  num2str(i)], 'fig');     
end
close all;
