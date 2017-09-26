function [ m2fit,m3fit] = grstat1( model1,model2,model3 )
%UNTITLED9 Summary of this function goes here
%   Detailed explanation goes here
%%% Compare spectrum of true signal and of estimated signal
%%%%%% NOTE highly dependent on noise used to generate model
% Define model inputs


    nelectrodes = size(model1.data,1); % number of electrodes
    f0 = model1.sampling_frequency;

 trials_true = model1.data;
 trials_spline = simulate_data(model2); 
 trials_standard = simulate_data(model3);

for electrode = 1:nelectrodes
    x=trials_true(electrode,:);
    x = x - ones(size(x)).*mean(x);
    TW = 2; % time-bandwidth product"
    [Sxx,~] = pmtm(x,TW,length(x),f0);
    

    %[~, Sxx] = mySpec(trials_true(electrode,:),f0,'noplot','tapers');
    h_true = Sxx;

      x=trials_standard(electrode,:);
    x = x - ones(size(x)).*mean(x);
    [Sxx, ~] = pmtm(x,TW,length(x),f0);
    %[~, Sxx] = mySpec(trials_standard(electrode,:),f0,'noplot','tapers');
    h_standard  = Sxx;

      x=trials_spline(electrode,:);
    x = x - ones(size(x)).*mean(x);
    [Sxx, ~] = pmtm(x,TW,length(x),f0);
   % [~, Sxx] = mySpec(trials_spline(electrode,:),f0,'noplot','tapers');
    h_spline  = Sxx;
   
    [H, X] = ecdf(h_true);
    [H1, X1] = ecdf(h_spline);
    [H2, X2] = ecdf(h_standard);
    
  

    ap = 2.2414; % for 95% confidence bounds
    total_observations = size(model1.data,2); %check! number of observations from which H is computed ? length of signal ?

    flag = 'biased'; % divide by 1/N

    R = xcov(trials_spline(electrode,:),flag); % autocovariance of estimated signal

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
    
    m2fit.xaxis = Q(1,:);
    m2fit.true = Q(2,:);
    m2fit.estimate = Q(4,:);
    m2fit.bound1 = Q(3,:);
    m2fit.bound2 = Q(5,:);
    
      if isnan(gr_spline)
        gr_stand=0;
    end
    m2fit.stat(electrode) = gr_spline;
%     figure;
%     plot(Q(1,:),Q(2,:),'k',Q(1,:),Q(4,:),'r',Q(1,:),Q(3,:),'--r',Q(1,:),Q(5,:),'--r')
%     hold on
    %%%%% standard calc
    R = xcov(trials_standard(electrode,:),flag); % autocovariance of estimated signal

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
    if isnan(gr_stand)
        gr_stand=0;
    end
    m3fit.xaxis = Q(1,:);
    m3fit.true = Q(2,:);
    m3fit.estimate = Q(4,:);
    m3fit.bound1 = Q(3,:);
    m3fit.bound2 = Q(5,:);
    m3fit.stat(electrode) = gr_stand;
  
end
% figure;
% plot(h_true,'k'); hold on; plot(h_standard,'g'); plot(h_spline,'r')


bfcorrection = 0.05/length(m2fit.stat);
m2fit.pvals =  pval2grstat(m2fit.stat,'grstatistic');
m3fit.pvals =  pval2grstat(m3fit.stat,'grstatistic');

m2fit.fails = find(m2fit.pvals<bfcorrection);
m3fit.fails = find(m3fit.pvals<bfcorrection);
end

