
ntrials = 2;
model_order_vals = [4 20 30 50];%[4 10 15 20 25 30 35 40 45 50 70 100];
jdist_results_spline = zeros(ntrials,length(model_order_vals));
time_results_spline = zeros(ntrials,length(model_order_vals));

jdist_results_standard = zeros(ntrials,length(model_order_vals));
time_results_standard = zeros(ntrials,length(model_order_vals));

accuracy_stand = zeros(ntrials,length(model_order_vals));
accuracy_spline = zeros(ntrials,length(model_order_vals));

config_spline;

%%% true adjacency matrix --------------------------------------
adj_true = model_true.true_coefficients;
adj_true(adj_true~=0)=1;
adj_true=sum(adj_true,3);
adj_true(adj_true~=0)=1;


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

        dj_stand = jdist(adj_true,adj_standard);
        dj_spline = jdist(adj_true,adj_mat);
        
        jdist_results_spline(ii,kk) = dj_spline;
        time_results_spline(ii,kk) = splinetime;
        
        jdist_results_standard(ii,kk) = dj_stand;
        time_results_standard(ii,kk) = standardtime;
        
        
        accuracy_stand(ii,kk) = network_accuracy(adj_true,adj_standard);
        accuracy_spline(ii,kk) = network_accuracy(adj_true,adj_mat);
    end
    fprintf([num2str(ii) ' \n']);
end



accspline = mean(accuracy_spline);
accstand = mean(accuracy_stand);

jd_stand = mean(jdist_results_standard);
jd_spline = mean(jdist_results_spline);
time_stand = mean(time_results_standard);
time_spline = mean(time_results_spline);

jd_stand_se = (1.96/sqrt(100)).*std(jdist_results_standard);
jd_spline_se = (1.96/sqrt(100)).*std(jdist_results_spline);
time_stand_se = (1.96/sqrt(100)).*std(time_results_standard);
time_spline_se = (1.96/sqrt(100)).*std(time_results_spline);
accspline_se = (1.96/sqrt(100)).*std(accuracy_spline);
accstand_se = (1.96/sqrt(100)).*std(accuracy_stand);

figure;
subplot 411
A=shadedErrorBar(model_order_vals,time_stand, time_stand_se , '-*b' , 1);
hold on;
B=shadedErrorBar(model_order_vals,time_spline, time_spline_se , '-*r' , 1);
h=legend([A.mainLine,B.mainLine],'Standard','Spline');
set(h,'FontSize',13)

title('Computation Time vs Lags','FontSize',15)
xlabel('Lag','FontSize',13)
ylabel('Computation Time (s)','FontSize',13)
% subplot 312
% A=shadedErrorBar(model_order_vals,1-jd_stand, jd_stand_se , '-*b' , 1);
% hold on;
% B=shadedErrorBar(model_order_vals,1-jd_spline, jd_spline_se , '-*r' , 1);
% title('Similarity vs Lags','FontSize',15)
% xlabel('Lag','FontSize',13)
% ylabel('Similarity','FontSize',13)
% h=legend([A.mainLine,B.mainLine],'Standard','Spline');
% set(h,'FontSize',13)


subplot 412
A=shadedErrorBar(model_order_vals,accstand, accstand_se , '-*b' , 1);
hold on;
B=shadedErrorBar(model_order_vals,accspline, accspline_se , '-*r' , 1);
title('Accuracy vs Lags','FontSize',15)
xlabel('Lag','FontSize',13)
ylabel('Accuracy','FontSize',13)
h=legend([A.mainLine,B.mainLine],'Standard','Spline');
set(h,'FontSize',13)



