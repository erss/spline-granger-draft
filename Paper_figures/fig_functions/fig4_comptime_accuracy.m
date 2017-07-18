%%%% lag v comp time & accuracy

%% Define loop parameters
ntrials = 100;
model_order_vals = [5:5:50];%[4 10 15 20 25 30 35 40 45 50 70 100];


for time = 1:3
          fprintf(['WINDOW SIZE: ' num2str(T(time)) ' \n']);

    T = [1 2 4];
    %% Declare results
    time_results_spline = zeros(ntrials,length(model_order_vals));
    time_results_standard = zeros(ntrials,length(model_order_vals));
    
    accuracy_stand = zeros(ntrials,length(model_order_vals));
    accuracy_spline = zeros(ntrials,length(model_order_vals));
    %% Set up sim
    % use nine node 20 lag in config : nine_node_order20_rdi;
    config_spline;
    model_true.T = T(time);   % time in seconds of window
    taxis = (1/model_true.sampling_frequency):(1/model_true.sampling_frequency):model_true.T;
    model_true.taxis = taxis;
    model_true.true_coefficients = nine_node_order20_rdi;
    %true adjacency matrix
    adj_true = model_true.true_coefficients;
    adj_true(adj_true~=0)=1;
    adj_true=sum(adj_true,3);
    adj_true(adj_true~=0)=1;
    
    %% Loop through every trial
    for ii = 1:ntrials
        
        %%% Simulate data
        model_true.data = simulate_data(model_true);
        
        
        for kk = 1:length(model_order_vals)
            
            model_order = model_order_vals(kk); % order used in model estimation
            model_true.estimated_model_order = model_order;
            if model_order == 4
                model_true.cntrl_pts = make_knots(4,3);
            elseif model_order == 10
                model_true.cntrl_pts  = make_knots(10,4);
            else
                model_true.cntrl_pts  = make_knots(model_order,floor(model_order/3));
            end
            
            tic
            [ adj_standard] = build_ar( model_true);
            standardtime  = toc;
            
            tic
            [ adj_mat] = build_ar_splines( model_true );
            splinetime  = toc;
            
            [accstand]= network_accuracy(adj_true,adj_standard);
            [accspline] = network_accuracy(adj_true,adj_mat);
            
            
            time_results_spline(ii,kk) = splinetime;
            time_results_standard(ii,kk) = standardtime;
            
            accuracy_stand(ii,kk) = accstand;
            accuracy_spline(ii,kk) = accspline;
        end
        fprintf([num2str(ii) ' \n']);
    end
    
    
    %% plot everything
    figure;
    subplot 211
    [mn, sem] = confidencebds(time_results_standard);
    A=shadedErrorBar(model_order_vals,mn, 2*sem, '-*b' , 1);
    hold on;
    [mn, sem] = confidencebds(time_results_spline);
    B=shadedErrorBar(model_order_vals,mn, 2*sem, '-*r' , 1);
    h=legend([A.mainLine,B.mainLine],'Standard','Spline');
    set(h,'FontSize',13)
    
    title('Computation Time vs Lags','FontSize',15)
    xlabel('Lag','FontSize',13)
    ylabel('Computation Time (s)','FontSize',13)
    
    subplot 212
    [mn, sem] = confidencebds(accuracy_stand);
    A=shadedErrorBar(model_order_vals,mn, 2*sem, '-*b' , 1);
    hold on;
    [mn, sem] = confidencebds(accuracy_spline);
    B=shadedErrorBar(model_order_vals,mn, 2*sem , '-*r' , 1);
    title('Accuracy vs Lags','FontSize',15)
    xlabel('Lag','FontSize',13)
    ylim([0 1])
    ylabel('Accuracy','FontSize',13)
    h=legend([A.mainLine,B.mainLine],'Standard','Spline');
    set(h,'FontSize',13)
    
    suptitle(num2str(T(time)));
end

%%% Save
h = get(0,'children');
for i=1:length(h)
    saveas(h(i), ['fig4_ctime'  num2str(i) 'lowfreq'], 'fig');
    
    
    
end
close all;
