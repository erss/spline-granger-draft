figure;
f0=model_true.sampling_frequency;
dt=1/f0;
T=model_true.T;

if strcmp(model_true.noise_type,'white')
    subplot (2,3,[1 3])
    plotchannels(dt:dt:T,model_true.data');
    
    subplot 234
    plotNetwork(model_true.network);
    
    subplot 235
    plotNetwork(model_standard.network);
    title(strcat({'Standard, '},num2str(model_standard.computation_time),{' s'},' Accuracy, ',num2str(model_standard.accuracy)))
    
    subplot 236
    plotNetwork(model_spline.network);
    title(strcat({'Spline, '},num2str(model_spline.computation_time),{' s'},' Accuracy, ',num2str(model_spline.accuracy)))
    
else
    
    subplot (2,2,[1 2])
    plotchannels(dt:dt:T,model_true.data');
    
    
    subplot 223
    plotNetwork(model_standard.network);
    title(strcat({'Standard, '},num2str(model_standard.computation_time),{' s'}))
    
    subplot 224
    plotNetwork(model_spline.network);
    title(strcat({'Spline, '},num2str(model_spline.computation_time),{' s'},' Overlap, ',num2str(model_spline.accuracy)))
end