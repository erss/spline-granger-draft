%%%% lag v comp time & accuracy

%% Define loop parameters
ntrials = 2;
model_order_vals = [1:3]; %[5:5:50];%[4 10 15 20 25 30 35 40 45 50 70 100];

%% Declare results

tp_spline = zeros(ntrials,length(model_order_vals));
tp_standard = zeros(ntrials,length(model_order_vals));


tn_spline = zeros(ntrials,length(model_order_vals));
tn_standard = zeros(ntrials,length(model_order_vals));

time_results_spline = zeros(ntrials,length(model_order_vals));
time_results_standard = zeros(ntrials,length(model_order_vals));

accuracy_stand = zeros(ntrials,length(model_order_vals));
accuracy_spline = zeros(ntrials,length(model_order_vals));
%% Set up sim
% use nine node 20 lag in config : nine_node_order20_rdi;
config_spline;
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

        [accstand, tpstand, tnstand, fnstand, fpstand]= network_accuracy(adj_true,adj_standard);
        [accspline, tpspline,tnspline,fnspline, fpspline] = network_accuracy(adj_true,adj_mat);
        
        
        time_results_spline(ii,kk) = splinetime;
        time_results_standard(ii,kk) = standardtime;
        
        tp_standard(ii,kk) = tpstand;
        tp_spline(ii,kk) = tpspline;
        
        tn_standard(ii,kk) = tnstand;
        tn_spline(ii,kk) = tnspline;
        
        
        fn_standard(ii,kk) = fnstand;
        fn_spline(ii,kk) = fnspline;
        
        fp_standard(ii,kk) = fpstand;
        fp_spline(ii,kk) = fpspline;
        
        
        
        accuracy_stand(ii,kk) = accstand;
        accuracy_spline(ii,kk) = accspline;
    end
    fprintf([num2str(ii) ' \n']);
end

%% calc means

mntpspline = mean(tp_spline,1); 
mntpstandard = mean(tp_standard,1); 

mnfpspline = mean(fp_spline,1); 
mnfpstandard = mean(fp_standard,1); 

mnfnspline = mean(fn_spline,1); 
mnfnstandard = mean(fn_standard,1); 

mntnspline = mean(tn_spline,1); 
mntnstandard =  mean(tn_standard,1); 

mnctspline = mean(time_results_spline,1);
mnctstandard = mean(time_results_standard,1); 

mnaccspline = mean(accuracy_spline,1);
mnaccstandard = mean(accuracy_stand,1);

%% calc se


semfpspline = std(fp_spline,1)/ntrials; 
semfpstandard = std(fp_standard,1)/ntrials; 
semfnspline = std(fn_spline,1)/ntrials; 
semfnstandard = std(fn_standard,1)/ntrials; 


semtpspline = std(tp_spline,1)/ntrials; 
semtpstandard = std(tp_standard,1)/ntrials; 

semtnspline = std(tn_spline,1)/ntrials; 
semtnstandard =  std(tn_standard,1)/ntrials; 

semctspline = std(time_results_spline,1)/ntrials;
semctstandard = std(time_results_standard,1)/ntrials; 

semaccspline = std(accuracy_spline,1)/ntrials;
semaccstandard = std(accuracy_stand,1)/ntrials;


%% plot everything
figure;
subplot 211
A=shadedErrorBar(model_order_vals,mnctstandard, 2*semctstandard, '-*b' , 1);
hold on;
B=shadedErrorBar(model_order_vals,mnctspline, 2*semctspline, '-*r' , 1);
h=legend([A.mainLine,B.mainLine],'Standard','Spline');
set(h,'FontSize',13)

title('Computation Time vs Lags','FontSize',15)
xlabel('Lag','FontSize',13)
ylabel('Computation Time (s)','FontSize',13)

subplot 212
A=shadedErrorBar(model_order_vals,mnaccstandard, 2*semaccstandard, '-*b' , 1);
hold on;
B=shadedErrorBar(model_order_vals,mnaccspline, 2*semaccspline , '-*r' , 1);
title('Accuracy vs Lags','FontSize',15)
xlabel('Lag','FontSize',13)
ylim([0 1])
ylabel('Accuracy','FontSize',13)
h=legend([A.mainLine,B.mainLine],'Standard','Spline');
set(h,'FontSize',13)

figure;
subplot 411
A=shadedErrorBar(model_order_vals,mntpstandard, 2*semtpstandard, '-*b' , 1);
hold on;
B=shadedErrorBar(model_order_vals,mntpspline, 2*semtpspline, '-*r' , 1);
title('TP vs Lags','FontSize',15)
xlabel('Lag','FontSize',13)
ylabel('TP','FontSize',13)
h=legend([A.mainLine,B.mainLine],'Standard','Spline');
set(h,'FontSize',13)


subplot 412
A=shadedErrorBar(model_order_vals,mntnstandard, 2*semtnstandard, '-*b' , 1);
hold on;
B=shadedErrorBar(model_order_vals,mntnspline, 2*semtnspline, '-*r' , 1);
title('TN vs Lags','FontSize',15)
xlabel('Lag','FontSize',13)
ylabel('TN','FontSize',13)
h=legend([A.mainLine,B.mainLine],'Standard','Spline');
set(h,'FontSize',13)


subplot 413
A=shadedErrorBar(model_order_vals,mnfpstandard, 2*semfpstandard, '-*b' , 1);
hold on;
B=shadedErrorBar(model_order_vals,mnfpspline, 2*semfpspline, '-*r' , 1);
title('FP vs Lags','FontSize',15)
xlabel('Lag','FontSize',13)
ylabel('FP','FontSize',13)
h=legend([A.mainLine,B.mainLine],'Standard','Spline');
set(h,'FontSize',13)

subplot 414
A=shadedErrorBar(model_order_vals,mnfnstandard, 2*semfnstandard, '-*b' , 1);
hold on;
B=shadedErrorBar(model_order_vals,mnfnspline, 2*semfnspline, '-*r' , 1);
title('FN vs Lags','FontSize',15)
xlabel('Lag','FontSize',13)
ylabel('FN','FontSize',13)
h=legend([A.mainLine,B.mainLine],'Standard','Spline');
set(h,'FontSize',13)

