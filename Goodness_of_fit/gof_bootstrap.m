function gof_bootstrap( model1,model2,model3)
%UNTITLED7 Summary of this function goes here
%   Detailed explanation goes here
% model1 = true
% model2 = spline
% model3 = standard
%%%% goodness_of_fit_bootstrap_coefficients

noise_type = model1.noise_type;


data = model1.data;
nelectrodes = size(data,1);
model_order = model1.estimated_model_order;

adj_true= model1.network;
adj_mat = model2.network;
if strcmp(noise_type,'white')
b = model1.true_coefficients;
nlags = size(b,3);  % true model order
end
bhat = model2.model_coefficients;
b_est_stand = model3.model_coefficients;

f0 = model1.sampling_frequency;


dt = 1/f0; % seconds
cntrl_pts = model1.cntrl_pts;

for electrode = 1:nelectrodes % plot fit for every electrode in network
 
    if sum(adj_mat(electrode,:))==0
       UB = zeros(model_order,nelectrodes);
       LB = zeros(model_order,nelectrodes);
    else
    [UB,LB]= myBootstrap(model2,electrode);
    end

%   figure;
  
      k=1;
    for i = 1:nelectrodes
        adj_sum = adj_mat + adj_true;
        adj_sum(adj_sum==2)=1;
        nconnections = sum(adj_sum(electrode,:));
        
        
   
        if (adj_mat(electrode,i) || adj_true(electrode,i))
       %     subplot(nconnections,1,k)
            plot(dt:dt:(model_order/f0),squeeze(real(bhat(electrode,i,:))),'r','LineWidth',1.5)
         % plot(squeeze(real(bhat(electrode,i,:))),'r','LineWidth',1.5)
               hold on
            plot(dt:dt:(model_order/f0),squeeze(real(b_est_stand(electrode,i,:))),'g','LineWidth',1.5)
          % plot(squeeze(real(b_est_stand(electrode,i,:))),'g','LineWidth',1.5)


             if strcmp(noise_type,'white') && nlags <= 5  
                 plot(dt:dt:(nlags/f0),squeeze(real(b(electrode,i,:))),'.k','MarkerSize',30);
                 %   plot(squeeze(real(b(electrode,i,:))),'.k','MarkerSize',30);
                % plot(dt:dt:(nlags/f0),squeeze(real(b(electrode,i,:))),'k','LineWidth',2);
             elseif strcmp(noise_type,'white')  %if nelectrodes > 1 && strcmp(noise_type,'white') 
                plot(dt:dt:(nlags/f0),squeeze(real(b(electrode,i,:))),'k','LineWidth',2);
            % plot(squeeze(real(b(electrode,i,:))),'k','LineWidth',2);

             end
             
             
            plot(cntrl_pts(2:end)./f0,squeeze(bhat(electrode,i,cntrl_pts(2:end))),'o')

            plot(dt:dt:(model_order/f0),real(LB(:,i)),'--r','LineWidth',1)
            plot(dt:dt:(model_order/f0),real(UB(:,i)),'--r','LineWidth',1)
% plot(squeeze(bhat(electrode,i,cntrl_pts(2:end))),'o')
% 
%             plot(real(LB(:,i)),'--r','LineWidth',1)
%             plot(real(UB(:,i)),'--r','LineWidth',1)
            
            
            
     %       str1 = strcat({'Influence of e'},num2str(i),{' on e'}, num2str(electrode));
     %       title(str1)
            xlabel('Lag (s)','FontSize',16);
            ylabel('Magnitude','FontSize',16);

            k= k +1;
            
        end

    end
    
    if ~strcmp(noise_type,'real')
        h = legend('Spline Estimated','Standard Estimated','True ');
        set(h,'FontSize',15,'Location','NorthEast');
    %    suptitle('Estimated Coefficient Fits');
    end
       if strcmp(noise_type,'real')
        h = legend('Spline Estimated','Standard Estimated');
        set(h,'FontSize',15,'Location','NorthEast');
       % suptitle('Estimated Coefficient Fits');
    end
    

end



end

