%%%% Get the data, estimate, and useful parameters.

data     = model_spline.data(1,model_spline.estimated_model_order+1:end);
estimate = model_spline.signal_estimate(1,:);
Fs       = model_spline.sampling_frequency;
N        = length(data);
T        = N/Fs;
m        = N/2;

t  = (0:N-1)/Fs;                    % time axis
fj = [(0:N/2-1)*(1/T)];             % freq axis, > 0

%%%% Compute the Fourier transform of estimate by hand.
X = zeros(1,length(fj));
n = (1:N);
for j=1:length(fj)
    X(j) = sum(estimate .* exp(-2*pi*1i*fj(j)*t));
end
%%%% Compute "model" power spectrum, scaled. See Priestley (6.1.24) and (6.2.4).
I_NX = 2/(N*4*pi) * (X.*conj(X));

%%%% Compute the Fourier transform of data by hand.
X = zeros(1,length(fj));
n = (1:N);
for j=1:length(fj)
    X(j) = sum(data .* exp(-2*pi*1i*fj(j)*t));
end
%%%% Compute "true" power spectrum, scaled. See Priestley (6.1.24) and (6.2.4).
f = 2/(N*4*pi) * (X.*conj(X));

%%%% Define frequency intervals of interest. See Priestley pg 487.
%%%% Here, use standard brain freq bands,
%     delta,theta,alpha,beta1, beta2, gamma,  NOTE: Ignore f > 50 Hz!!!!
%Ij = [0,4; 4,8;  8,12; 12,20; 20,30; 30,50]; %30,80; 80,120]; %120,fj(end)];
Ij = [0,4; 4,8;  8,12; 12,20; 20,30; 30,50];
k  = size(Ij,1);
Sj = zeros(k,1);
nj = zeros(k,1);

%%%% Compute the summed, scaled periodogram ordinates for each interval. See Priestley (6.2.179).
for j=1:k
    indices_for_f_in_interval = find(fj>Ij(j,1) & fj<=Ij(j,2));
    nj(j) = length(indices_for_f_in_interval);
    Sj(j) = 1/(2*nj(j)) * sum( I_NX(indices_for_f_in_interval) ./ f(indices_for_f_in_interval));
end

%%%% Compute the test statistic. See Priestley (6.2.180).
C = 1 + 1/(3*(k-1)) * ( sum(1/(2*nj)) - 1/(2*sum(nj)));
M = C * ( 2*sum(nj) * log( sum(nj.*Sj)/sum(nj) ) - 2*sum( nj.*log(Sj) ) );

%%%% Under the null hypothesis that the specified form of f is the true
%%%% one. Compute chi2 with k-1 d.o.f.
p=1-chi2cdf(M,k-1);

fprintf(['p value is ' num2str(p,3) '\n'])

%%%% Useful plots.

%%%% Plot the data, and the estimate.
subplot(2,1,1)
plot(t,data, 'k')
hold on
plot(t,estimate, 'r')
hold off
axis tight
xlabel('Time [s]')

subplot(2,2,3)
plot(fj, f, 'k')
hold on
plot(fj, I_NX, 'r')
hold off
xlabel('Frequency [Hz]')

subplot(2,2,4)
x = (0:0.1:40);
chi2 = zeros(size(x));
for j=1:length(x)
    chi2(j) = chi2pdf(x(j),k-1);
end
plot(x,chi2)
hold on
plot([M,M], [0, max(chi2)], 'r')
hold off
xlabel('M values')
ylabel('chi2')
title(['p value is ' num2str(p,3)])


