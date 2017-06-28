figure;
f0=model_true.sampling_frequency;
dt=1/f0;
T=model_true.T;
subplot (2,3,[1 3])
plotchannels(dt:dt:T,model_true.data');

subplot 234
plotNetwork(model_true.network);


subplot 235
plotNetwork(model_standard.network);
title(strcat({'Standard, '},num2str(model_standard.computation_time),{' s'},' Similarity, ',num2str(model_standard.jaccard_similarity)))

subplot 236
plotNetwork(model_spline.network);
title(strcat({'Spline, '},num2str(model_spline.computation_time),{' s'},' Similarity, ',num2str(model_spline.jaccard_similarity)))

