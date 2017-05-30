%%%% goodness_of_fit_bootstrap_coefficients
if ~exist('noise_type','var')
  noise_type = 'white';
end

for electrode = 1:nelectrodes % plot fit for every electrode in network
 
    if sum(adj_mat(electrode,:))==0
       UB = zeros(model_order,nelectrodes);
       LB = zeros(model_order,nelectrodes);
    else
    [UB,LB]= myBootstrap(data,adj_mat,model_order,electrode,cntrl_pts);
    end;

    figure;
  
      k=1;
    for i = 1:nelectrodes
        adj_sum = adj_mat + adj_true;
        adj_sum(adj_sum==2)=1;
        nconnections = sum(adj_sum(electrode,:));
        
        
   
        if (adj_mat(electrode,i) || adj_true(electrode,i))
            subplot(1,nconnections,k)
            plot(dt:dt:(model_order/f0),squeeze(real(bhat(electrode,i,:))),'r','LineWidth',1.5)
               hold on


             if nelectrodes == 1 && strcmp(noise_type,'white')
                plot(dt:dt:(nlags/f0),squeeze(real(b(electrode,i,:))),'.k','MarkerSize',30);
             elseif nelectrodes > 1 && strcmp(noise_type,'white')
                plot(dt:dt:(nlags/f0),squeeze(real(b(electrode,i,:))),'k','LineWidth',1.5);
             end
             
             
            plot(cntrl_pts(2:end)./f0,squeeze(bhat(electrode,i,cntrl_pts(2:end))),'o')

            plot(dt:dt:(model_order/f0),real(LB(:,i)),'--r','LineWidth',1)
            plot(dt:dt:(model_order/f0),real(UB(:,i)),'--r','LineWidth',1)

            
            
            
            str1 = strcat({'Influence of e'},num2str(i),{' on e'}, num2str(electrode));
            title(str1)
            xlabel('Lag (s)','FontSize',14);
            ylabel('Magnitude','FontSize',14);

            k= k +1;
            
        end

    end


    h = legend('Spline Estimated Coefficients','True AR coefficients');
    set(h,'FontSize',14,'Location','SouthEast');
    suptitle('Estimated Coefficient Fits');


end

