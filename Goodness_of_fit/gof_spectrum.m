function [ output_args ] = gof_spectrum( model1,model2,model3 )
%UNTITLED9 Summary of this function goes here
%   Detailed explanation goes here
%%% Compare spectrum of true signal and of estimated signal
%%%%%% NOTE highly dependent on noise used to generate model
% Define model inputs

nrealizations = model1.nrealizations;
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

    lb = round(0.025*nrealizations);
    if lb == 0
       lb = 1; 
    end
    ub = round(0.975*nrealizations);
    
    figure;
    plot(X,H,'k','LineWidth',2);
    hold on;
    plot(X1,H1,'r','LineWidth',2)
    plot(X2,H2,'g','LineWidth',2)
    
%    plot(X,true_H0(lb,:),'--k','LineWidth',2);
    plot(X1,spline_H0(lb,:),'--r','LineWidth',2)
    plot(X2,standard_H0(lb,:),'--g','LineWidth',2)
%    plot(X,true_H0(ub,:),'--k','LineWidth',2);
    plot(X1,spline_H0(ub,:),'--r','LineWidth',2)
    plot(X2,standard_H0(ub,:),'--g','LineWidth',2)
    
end

end
