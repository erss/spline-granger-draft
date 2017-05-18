%%%% goodness_of_fit_bootstrap_coefficients

for electrode = 1:nelectrodes % plot fit for every electrode in network
    nsurrogates = 100; % number of surrogates

    [b_surrogates, bounds]= myBootstrap(data,adj_mat,model_order,electrode,nsurrogates,cntrl_pts);


    figure;
      j=1; 
    for i = 1:nelectrodes
        subplot(1,nelectrodes,i)
       
    
        plot(dt:dt:(model_order/f0),squeeze(real(bhat(electrode,i,:))),'r','LineWidth',1.5)
           hold on
        plot(cntrl_pts(2:end)./f0,squeeze(bhat(electrode,i,cntrl_pts(2:end))),'o')
        plot(dt:dt:(model_order/f0),squeeze(real(bounds(1,j:j+model_order-1))),'--r')
        plot(dt:dt:(model_order/f0),squeeze(real(bounds(2,j:j+model_order-1))),'--r')
        
         if nelectrodes == 1
            plot(dt:dt:(nlags/f0),squeeze(real(b(electrode,i,:))),'.k','MarkerSize',30);
        else
            plot(dt:dt:(nlags/f0),squeeze(real(b(electrode,i,:))),'k','LineWidth',1.5);
        end
            str1 = strcat({'Influence of e'},num2str(i),{' on e'}, num2str(electrode));
        title(str1)
    xlabel('Lag (s)','FontSize',14);
    ylabel('Magnitude','FontSize',14);
    j= j+model_order;
    end


    h = legend('Real AR Coefficients','Spline Estimated Coefficients');
    set(h,'FontSize',14,'Location','SouthEast');
    suptitle('Estimated Coefficient Fits');


end

