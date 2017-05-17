%%%% goodness_of_fit_bootstrap_coefficients

%electrode = 1; % which electrode in network to generate data for 

for electrode = 1:nelectrodes
    nsurrogates = 100; % number of surrogates

    [b_surrogates, bounds]= myBootstrap(data,adj_mat,model_order,electrode,nsurrogates,cntrl_pts);


%     bounds_plus = zeros(1,nelectrodes,model_order);
%     bounds_minus = zeros(1,nelectrodes,model_order);
% 
%     j=1;   
%     for p = 1:nelectrodes
%         if adj_mat(electrode,p) == 1
%             bounds_plus(electrode,p,:) = bounds(1,j:j+model_order-1);
%             bounds_minus(electrode,p,:) = bounds(2,j:j+model_order-1);
%             j= j+model_order;
%         end    
%     end

    figure;
      j=1; 
    for i = 1:nelectrodes
        subplot(1,nelectrodes,i)
        plot(squeeze(real(b(electrode,i,:))),'*k');
        hold on
        plot(squeeze(real(bounds(1,j:j+model_order-1))),'--r')
        plot(squeeze(real(bounds(2,j:j+model_order-1))),'--r')

%         plot(squeeze(real(bounds_plus(electrode,i,:))),'--r');
%         plot(squeeze(real(bounds_minus(electrode,i,:))),'--r');
            xlabel('Lag','FontSize',14);
    ylabel('Magnitude','FontSize',14);
    j= j+model_order;
    end


    h = legend('Real AR Coefficients','Spline Estimated Coefficients');
    set(h,'FontSize',14,'Location','SouthEast');
    suptitle('Estimated Coefficient Fits');


end

