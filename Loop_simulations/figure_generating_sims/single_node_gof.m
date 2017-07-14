%%%%% single node analysis
%%% set model coefficitens to single_node_order20 & single_node_low_freq
   clear all;
for kk = 2
 
    ntrials = 50;
    fails =[];
    for i = 1:ntrials
        config_spline;
        if kk==1
            model_true.true_coefficients = single_node_low_freq; %%%% MODIFY COEFFICIENTS HERE!
        else
            model_true.true_coefficients = single_node_order20; %%%% MODIFY COEFFICIENTS HERE!
            
        end
        model_true.model_coefficients = model_true.true_coefficients;
        simulate_network;
        infer_network;
        
        ct_spline(i)  = model_spline.computation_time;
        ct_standard(i) = model_standard.computation_time;
        
        
        acc_spline(i)  = model_spline.accuracy;
        acc_standard(i) = model_standard.accuracy;
        
        [nw dwstandard(i)] = gof_residuals(model_standard);
        fails = [fails nw];
        [nw dwspline(i)] = gof_residuals(model_spline);
        fails = [fails nw];
        
        [ts_spline(i), ts_stand(i)] = gof_spectrum(model_true,model_spline,model_standard);
        close all;
        
        fprintf(num2str(i));
    end
    %%
    figure;
    subplot(3,2,1)
    plot(model_true.taxis,model_true.data(1,:),'k');
    title('Simulated Signal','FontSize',14);
    xlabel('Time (s)','FontSize',14);
    subplot(3,2,2)
    mySpec(model_true.data(1,:),model_true.sampling_frequency,'yesplot','tapers');
    subplot(3,2,[3 4])
    gof_bootstrap(model_true,model_spline,model_standard);
    title('AR Coefficients','FontSize',14);
    subplot(3,2,5)
    gof_spectrum(model_true,model_spline,model_standard);
    title('Integrated Spectrum Test','FontSize',14);
    subplot(3,2,6)
    data  = model_true.data;
    
    nelectrodes = size(data,1);
    model_order = model_true.estimated_model_order;
    datap = data(:,model_order+1:end);
    yestimate = model_spline.signal_estimate;
    
    residuals = datap-yestimate;
    autocorr(residuals);
    title('Autocorrelation of Residuals','FontSize',14);
    xlabel('Lag Index','FontSize',14);
    ylabel('Autocorrelation','FontSize',14);
    
    %%

    
     Labels = {'Standard', 'Spline'};
    figure;
    subplot 131
barplot(Labels,ct_standard,ct_spline)
    ylabel('Computation time (s)','FontSize',14)
    
    subplot 132
    barplot(Labels,ts_stand,ts_spline)

     xlim=get(gca,'xlim');
    hold on
    plot(xlim,[2.2414 2.2414],'--r','LineWidth',2)
    ylabel('Grenander & Rosenblatt Statistic','FontSize',14)
    
    subplot 133
    barplot(Labels,dwstandard,dwspline)
    xlim=get(gca,'xlim');
    hold on
    plot(xlim,[1 1],'--r','LineWidth',2)
    plot(xlim,[3 3],'--r','LineWidth',2)
    ylabel('Durbin-Watson Statistic','FontSize',14)
    
%     h = get(0,'children');
%     for i=1:length(h)
%         if kk == 1
%             saveas(h(i), ['single_node_gof'  num2str(i) 'lowfreq'], 'fig');
%         else
%             saveas(h(i), ['single_node_gof'  num2str(i) 'specpeak'], 'fig');
%         end
%         
%     end
%     close all;
    
    
end

