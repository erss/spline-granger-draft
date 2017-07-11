function [ cons ] = gof_consistency( model_spline)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

ntrials = 1;
order = model_spline.estimated_model_order;
data = model_spline.data(:,order+1:end);
nelectrodes = size(data,1);
N = size(data,2);
trials_true     = zeros(nelectrodes,N,ntrials);
trials_spline     = zeros(nelectrodes,N,ntrials);


for i = 1:ntrials
    dataspline = simulate_data(model_spline);
    trials_spline(:,:,i) =  model_spline.signal_estimate; %dataspline(order+1:end);
    trials_true(:,:,i) = data;
end
s = ntrials*size(data,2);
Y=reshape(trials_spline,[nelectrodes s]); % stack trials
X=reshape(trials_true,[nelectrodes s]);

Rr = (X*X')/(s-1);           % covariance estimate of real data
Rs = (Y*Y')/(s-1);           % covariance estimate of estimated spline data

cons = 1 - norm(Rs-Rr)/norm(Rr); % compare matrix norms


end

