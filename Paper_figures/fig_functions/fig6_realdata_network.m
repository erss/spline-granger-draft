%%%% Real Data!
config_spline;
model_true.noise_type = 'real'; % 'white', 'pink', 'real'
 model_true.sztype = 'presz'; % presz
% model_true.ntwk = [2 7 18 22 42 46 90 80 77 ];     % badchannels = [1,9,21,32,83, 8,31];

 badchannels = [1,8,9,13,20,21,22,23,24,25,26,31,32,34,38,41,48,49,50,68,69,71,77,78,79,82,83,87,88 ,89];
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

h = get(0,'children');
for i=1:length(h)
saveas(h(i), ['fig6_ntwk'  num2str(i)], 'fig');     
end
close all;