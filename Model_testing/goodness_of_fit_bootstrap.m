%%%% goodness_of_fit_bootstrap_coefficients

for electrode = 1:nelectrodes % plot fit for every electrode in network
    nsurrogates = 100; % number of surrogates

    [b_surrogates, bounds]= myBootstrap(data,adj_mat,model_order,electrode,nsurrogates,cntrl_pts);


    figure;
      j=1; 
    for i = 1:nelectrodes
        subplot(1,nelectrodes,i)
        plot(squeeze(real(b(electrode,i,:))),'*k');
        hold on
        plot(squeeze(real(bounds(1,j:j+model_order-1))),'--r')
        plot(squeeze(real(bounds(2,j:j+model_order-1))),'--r')

    xlabel('Lag','FontSize',14);
    ylabel('Magnitude','FontSize',14);
    j= j+model_order;
    end


    h = legend('Real AR Coefficients','Spline Estimated Coefficients');
    set(h,'FontSize',14,'Location','SouthEast');
    suptitle('Estimated Coefficient Fits');


end

