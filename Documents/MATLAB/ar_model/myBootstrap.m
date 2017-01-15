function [ conf_int ] = myBootstrap( data, nshuffles, nsamples )
% myBootstrap takes DATA vector and draws NSAMPLES number of values
% randomly with replacement. For each sample, it computes the sample mean.
% Repeats this NSHUFFLES number of times and from this computes 95%
% confidence intervals.


 for ii = 1:nshuffles  % shuffles data nshuffles number of times
    dsample = datasample(data,nsamples);  % randomly samples nsamples with replacement
    meansample(ii) = mean(dsample);       % takes mean for every sample
 end
 
 std_means = std(meansample); % standard deviation of means
 avg_means = mean(meansample); % mean of means
 
 
 meansample = -1.*meansample;
 meansample = sort(meansample);   % sort data
 lower_ind = nshuffles * 0.025;   % finds index of 2.5%
 upper_ind = nshuffles * 0.975;   %            ... 97.5 %
 
 
 %histogram(meansample)
 % if lower_ind and upper_ind aren't integer numbers, take average of
 % ceiling and floor index of mean sample
 CI(1) = mean([meansample(ceil(lower_ind)),meansample(floor(lower_ind))]);
 CI(2) = mean([meansample(ceil(upper_ind)),meansample(floor(upper_ind))]);
 
 % 95% confidence interval
 conf_int = [avg_means - std_means*(CI(1)-CI(2)); avg_means + std_means*(CI(1)-CI(2))];
 

end

