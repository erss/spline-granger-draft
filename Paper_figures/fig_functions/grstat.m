function [ m2fit,m3fit] = grstat( model1,model2,model3 )
%UNTITLED9 Summary of this function goes here
%   Detailed explanation goes here
%%% Compare spectrum of true signal and of estimated signal
%%%%%% NOTE highly dependent on noise used to generate model
% Define model inputs

nrealizations = 100;%model1.nrealizations;
noise_type = model1.noise_type;

    nelectrodes = size(model1.data,1); % number of electrodes
    T = model1.T;      % total length of recording (seconds)
    f0 = model1.sampling_frequency;
    N=T*f0;


trials_true     = zeros(nelectrodes,N,nrealizations);
trials_spline   =  zeros(nelectrodes,N,nrealizations);
trials_standard =  zeros(nelectrodes,N,nrealizations);

for i = 1:nrealizations 
    
    if strcmp(noise_type,'white')
        trials_true(:,:,i) = simulate_data(model1);
    end
        trials_spline(:,:,i) = simulate_data(model2); 
        trials_standard(:,:,i) = simulate_data(model3);
       
end

if ~strcmp(noise_type,'white')
    trials_true = model1.data;
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
    mat2 = squeeze(trials_standard(electrode,:,:))';
    mat3 = squeeze(trials_spline(electrode,:,:))';
    
    
    standard_spectrum = [];
    spline_spectrum = [];
    
    
    standard_H = [];
    spline_H = [];
    for i = 1:nrealizations
        
        if strcmp(noise_type,'white')
            [~, Sxx] = mySpec(mat1(i,:),f0,'noplot','notapers');
            true_spectrum = [true_spectrum; Sxx];
            true_H = [true_H; ecdf(Sxx)'];
        end
        
        [~, Sxx] = mySpec(mat2(i,:),f0,'noplot','notapers');
        standard_spectrum = [standard_spectrum; Sxx];
        standard_H = [standard_H; ecdf(Sxx)'];
        
        [~, Sxx] = mySpec(mat3(i,:),f0,'noplot','notapers');
        spline_spectrum = [spline_spectrum; Sxx];
        spline_H = [spline_H; ecdf(Sxx)'];
    end
    
    
    h_standard = mean(standard_spectrum);
    h_spline = mean(spline_spectrum);
     if strcmp(noise_type,'white')
         h_true = mean(true_spectrum);
         [H, X] = ecdf(h_true);
         true_H0 = sort(true_H,1);
     end
    [H1, X1] = ecdf(h_spline);
    [H2, X2] = ecdf(h_standard);
    
    
    spline_H0 = sort(spline_H,1);
    standard_H0 = sort(standard_H,1);


    ap = 2.2414; % for 95% confidence bounds
    total_observations = size(model1.data,2); %check! number of observations from which H is computed ? length of signal ?

    flag = 'biased'; % divide by 1/N

    R = xcov(trials_spline(:,:,1),flag); % autocovariance of estimated signal

    G = sum(R(970:end-30).^2);
    G = G/(4*pi);

    conf1 = ap*sqrt(8*pi*G/total_observations);
    UB = H1 + conf1;
    LB = H1 - conf1;
     %%% spline calc
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
    gr_spline = max(sqrt(total_observations)*abs(Qp(2,:)-Qp(4,:)));
%     figure;
%     plot(Q(1,:),Q(2,:),'k',Q(1,:),Q(4,:),'r',Q(1,:),Q(3,:),'--r',Q(1,:),Q(5,:),'--r')
%     hold on
   m2fit.xaxis = Q(1,:);
    m2fit.true = Q(2,:);
    m2fit.estimate = Q(4,:);
    m2fit.bound1 = Q(3,:);
    m2fit.bound2 = Q(5,:);
    m2fit.stat = gr_spline;
    %%%%% standard calc
    R = xcov(trials_standard(:,:,1),flag); % autocovariance of estimated signal

    G = sum(sum(R(970:end)).^2);
    G = G/(4*pi);

    conf2 = ap*sqrt(8*pi*G/total_observations);
    UB = H2 + conf2;
    LB = H2 - conf2;
       
    q1 = [X' X2'];
    Q = NaN(5,length(q1));
    Q(1,:) = q1;
    Q(2,:) = interp1q(X,H,Q(1,:)')';
    Q(3,:) = interp1q(X2,LB,Q(1,:)')';
    Q(4,:) = interp1q(X2,H1,Q(1,:)')';
    Q(5,:) = interp1q(X2,UB,Q(1,:)')'';

    Q = sortrows(Q',1)';
    Qp=Q';
    Qp = Qp(~any(isnan(Qp),2),:)';
    
   % plot(Q(1,:),Q(4,:),'g',Q(1,:),Q(3,:),'--g',Q(1,:),Q(5,:),'--g')
    % pg 476
    gr_stand = max(sqrt(total_observations)*abs(Qp(2,:)-Qp(4,:)));
       m3fit.xaxis = Q(1,:);
    m3fit.true = Q(2,:);
    m3fit.estimate = Q(4,:);
    m3fit.bound1 = Q(3,:);
    m3fit.bound2 = Q(5,:);
    m3fit.stat = gr_stand;
end
% figure;
% plot(h_true,'k'); hold on; plot(h_standard,'g'); plot(h_spline,'r')
end

