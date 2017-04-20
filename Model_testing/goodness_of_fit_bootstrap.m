%%%% goodness_of_fit_bootstrap_coefficients

electrode = 1; % which electrode in network to generate data for 
nsurrogates = 100; % number of surrogates

[b_surrogates, bounds]= myBootstrap(data,adj_mat,model_order,electrode,nsurrogates,cntrl_pts);


bounds_plus = zeros(1,nelectrodes,model_order);
bounds_minus = zeros(1,nelectrodes,model_order);

j=1;   
for p = 1:nelectrodes
    if adj_mat(electrode,p) == 1
        bounds_plus(electrode,p,:) = bounds(1,j:j+model_order-1);
        bounds_minus(electrode,p,:) = bounds(2,j:j+model_order-1);
        j= j+model_order;
    end    
end

for i = 1:nelectrodes
    figure;
    plot(squeeze(real(b(electrode,i,:))),'k');
    hold on
    plot(squeeze(real(bounds_plus(electrode,i,:))),'--r');
    plot(squeeze(real(bounds_minus(electrode,i,:))),'--r');
    
end
      



