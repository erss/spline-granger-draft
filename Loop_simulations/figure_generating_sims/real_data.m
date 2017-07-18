%%%% Real Data!
for i = 1%:90
config_spline;

 model_true.sztype = 'sz'; % presz
 model_true.ntwk = 22 %[2 7 18 22 42 46 90 80 77 ]; % i;%[77];     % badchannels = [1,9,21,32,83, 8,31];
%  [2 7 18 22 42 46 90 80 77 ]; 
 simulate_network;
 infer_network;

 
 figure;
subplot 221

plot(model_true.taxis,model_true.data,'k','LineWidth',2)
     hold on;
     plot([0,.2], [min(model_true.data(1,:)), min(model_true.data(1,:))], 'r', 'LineWidth', 2);
set(gca,'YTickLabel',[]);
set(gca,'XTickLabel',[]);
set(gca,'XTick',[]);
set(gca,'YTick',[]);

ylabel('Signal','FontSize',18);
xlabel('Time (s)','FontSize',18);

title('Single Electrode Recording','FontSize',20);
%title(num2str(i),'FontSize',20);

subplot 222

list=lines(size(model_true.data,1));
for i =1:size(model_true.data,1)
    mySpec(model_true.data(i,:),model_true.sampling_frequency,'yesplot','notapers','k');
    hold on;
  
end
title('Spectrogram','FontSize',20);

% 
 subplot(2,2,[3 4])
 gof_bootstrap(model_true,model_spline,model_standard);
title('Estimated Coefficient Fits','FontSize',20);
end

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

