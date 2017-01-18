function [ conf_int ] = myBootstrap( data, b, nshuffles, nsamples )
% myBootstrap takes DATA vector and draws NSAMPLES number of values
% randomly with replacement. For each sample, it computes the sample mean.
% Repeats this NSHUFFLES number of times and from this computes 95%
% confidence intervals.


 for ii = 1:nshuffles  % shuffles data nshuffles number of times
    dsample = datasample(data,nsamples);  % randomly samples nsamples with replacement
    delta(ii) = mean(dsample) - b;       % takes mean for every sample and 
                                        % calculates difference from 'true'
 end
 
 delta = sort(delta);   % sort data
 lower_ind = nshuffles * 0.025;   % finds index of 2.5%
 upper_ind = nshuffles * 0.975;   %            ... 97.5 %
 
 
 %histogram(delta)
 % if lower_ind and upper_ind aren't integer numbers, take average
 % meansample at ceiling and floor index
 %CI(1) = mean([meansample(ceil(lower_ind)),meansample(floor(lower_ind))]);
%  %CI(2) = mean([meansample(ceil(upper_ind)),meansample(floor(upper_ind))]);
%  CI(1) = meansample(lower_ind);
%  CI(2) = meansample(upper_ind);

 % 95% confidence interval
 conf_int = [b- delta(upper_ind); b - delta(lower_ind)];
 

end

