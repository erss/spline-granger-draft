% Bartlett homogeneity of chi squared test.
%
%NOTES.
% x. OPTION - Do not use "Fourier frequencies". Instead use 1 Hz freq resolution?

function [result] = gof_bartlett(model, Ij)

    n_sensors = size(model.data,1);
    result    = cell(n_sensors,1);
    
    for nk=1:n_sensors

        data     = model.data(nk,model.estimated_model_order+1:end);
        estimate = model.signal_estimate(nk,:);
        Fs       = model.sampling_frequency;
        N        = length(data);
        T        = N/Fs;

        t  = (0:N-1)/Fs;                    % time axis
        fj = [(0:N/2-1)*(1/T)];            % freq axis, > 0
        %fj = (0:1:(N/2-1)/T);               % Subsample, every 1 Hz.

        %%%% Compute the Fourier transform of estimate by hand.
        X = zeros(1,length(fj));
        for j=1:length(fj)
            X(j) = sum(estimate .* exp(-2*pi*1i*fj(j)*t));
        end
        %%%% Compute "model" power spectrum, scaled. See Priestley (6.1.24) and (6.2.4).
        I_NX = 2/(N*4*pi) * (X.*conj(X));

        %%%% Compute the Fourier transform of data by hand.
        X = zeros(1,length(fj));
        for j=1:length(fj)
            X(j) = sum(data .* exp(-2*pi*1i*fj(j)*t));
        end
        %%%% Compute "true" power spectrum, scaled. See Priestley (6.1.24) and (6.2.4).
        f = 2/(N*4*pi) * (X.*conj(X));

        %%%% Define frequency intervals of interest. See Priestley pg 487.
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
        
        %%%% Save the results.
        result{nk}.p = p;
        result{nk}.M = M;
        result{nk}.k = k;
        result{nk}.t = t;
        result{nk}.data = data;
        result{nk}.estimate = estimate;
        result{nk}.fj = fj;
        result{nk}.f  = f;
        result{nk}.I_NX = I_NX;
    end
end