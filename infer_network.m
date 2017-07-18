%%% True network ---------------------------------------------------------

if strcmp(model_true.noise_type,'white')
    adj_true = model_true.true_coefficients;
    adj_true(adj_true~=0)=1;
    adj_true=sum(adj_true,3);
    adj_true(adj_true~=0)=1;
    
    
  %   adj_true(1:size(adj_true,1)+1:end) =0;
    
    model_true.network = adj_true;
    

else
    nelectrodes = size(model_true.data,1);
    adj_true = ones(nelectrodes,nelectrodes);
    model_true.network = adj_true;
end

%%% Fit standard AR to data ----------------------------------------------

tic
[ adj_standard] = build_ar( model_true );
standardtime  = toc;
[ bhat, yhat] = estimate_standard( model_true, adj_standard);


model_standard = model_true;
model_standard.model_coefficients = bhat;
model_standard.computation_time = standardtime;
model_standard.signal_estimate = yhat;
model_standard.network = adj_standard;

if strcmp(model_true.noise_type,'white')
model_standard.jaccard_similarity = 1- jdist(model_true.network,model_standard.network);
model_standard.accuracy = network_accuracy(model_true.network,model_standard.network);
end

%%% Fit spline to data ---------------------------------------------------
tic
[ adj_spline] = build_ar_splines( model_true);
splinetime  = toc;
[ bhat, yhat] = estimate_coefficient_fits( model_true, adj_spline);


model_spline = model_true;
model_spline.model_coefficients = bhat;
model_spline.computation_time = splinetime;
model_spline.signal_estimate = yhat;
model_spline.network = adj_spline;

if strcmp(model_true.noise_type,'white')
model_spline.jaccard_similarity = 1- jdist(model_true.network,model_spline.network);
model_spline.accuracy = network_accuracy(model_true.network,model_spline.network);
end

if ~strcmp(model_true.noise_type,'white')
    model_spline.jaccard_similarity = 1- jdist(model_standard.network,model_spline.network);
    model_spline.accuracy = network_accuracy(model_standard.network,model_spline.network);
    
    model_standard.jaccard_similarity = 1- jdist(model_standard.network,model_spline.network);
    model_standard.accuracy = network_accuracy(model_standard.network,model_spline.network);
end
% model_spline.covb = covariance_b;
% model_spline.design_matrix = dm;
