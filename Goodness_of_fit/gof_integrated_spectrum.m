function [ gr_estimated] = gof_integrated_spectrum( model_true,model_estimated )
% GOF_INTEGRATED_SPECTRUM computes 

nrealizations = model_true.nrealizations;
noise_type = model_true.noise_type;

    nelectrodes = size(model_true.data,1); % number of electrodes
    T = model_true.T;      % total length of recording (seconds)
    f0 = model_true.sampling_frequency;
    N=T*f0;


trials_true     = zeros(nelectrodes,N,nrealizations);
trials_estimated   =  zeros(nelectrodes,N,nrealizations);

for i = 1:nrealizations 
    
    if strcmp(noise_type,'white')
        trials_true(:,:,i) = simulate_data(model_true);
    end
        trials_estimated(:,:,i) = simulate_data(model_estimated);        
end

if ~strcmp(noise_type,'white')
    trials_true = model_true.data;
end


for electrode = 1:nelectrodes
    
    if strcmp(noise_type,'white')
        mat1 = squeeze(trials_true(electrode,:,:))';
        true_spectrum = [];
        true_H = [];
    else
            [~, Sxx] = mySpec(trials_true(electrode,:),f0,'noplot','notapers');
            true_spectrum = Sxx;
            [H, X] = ecdf(true_spectrum);
            H = H';
        
    end
    mat3 = squeeze(trials_estimated(electrode,:,:))';
    
    
    estimated_spectrum = [];
    estimated_cdf = [];

    for i = 1:nrealizations
        
        if strcmp(noise_type,'white')
            [~, Sxx] = mySpec(mat1(i,:),f0,'noplot','notapers');
            true_spectrum = [true_spectrum; Sxx];
            true_H = [true_H; ecdf(Sxx)'];
        end
        
        [~, Sxx] = mySpec(mat3(i,:),f0,'noplot','notapers');
        estimated_spectrum = [estimated_spectrum; Sxx];
        estimated_cdf = [estimated_cdf; ecdf(Sxx)'];
        
    end
    
    
    h_estimated = mean(estimated_spectrum);
     if strcmp(noise_type,'white')
         h_true = mean(true_spectrum);
         [H, X] = ecdf(h_true);
     end
    [H1, X1] = ecdf(h_estimated);
    

   plot(X,H,'k','LineWidth',2);
   hold on;
   plot(X1,H1,'--r','LineWidth',2)
    

    ap = 2.2414; % for 95% confidence bounds
    total_observations = size(model_true.data,2); %check! number of observations from which H is computed ? length of signal ?

    flag = 'biased'; % divide by 1/N

    R = xcov(trials_estimated(:,:,1),flag); % autocovariance of estimated signal

    G = sum(R(970:end-30).^2);
    G = G/(4*pi);

    conf1 = ap*sqrt(8*pi*G/total_observations);
    UB = H1 + conf1;
    LB = H1 - conf1;
    
  
    plot(X1,UB, '--r');
    plot(X1,LB, '--r');
    axis tight
    
    
    h=legend('True Network','Estimated Network');
    set(h,'FontSize',15,'Location','SouthEast');
    %%% Compute amount of time in confidence bounds

    q1 = [X' X1'];
    Q = NaN(5,length(q1));
    Q(1,:) = q1;
    Q(2,:) = interp1q(X,H,Q(1,:)')';
    Q(3,:) = interp1q(X1,LB,Q(1,:)')';
    Q(4,:) = interp1q(X1,H1,Q(1,:)')';
    Q(5,:) = interp1q(X1,UB,Q(1,:)')'';
    
    Q = sortrows(Q',1)';
    Qp=Q';
    Qp = Qp(~any(isnan(Qp),2),:)';
    % pg 476
    gr_estimated = max(sqrt(total_observations)*abs(Qp(2,:)-Qp(4,:)));
    
 
end

end

