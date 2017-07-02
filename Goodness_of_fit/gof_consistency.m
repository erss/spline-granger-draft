function [ cons ] = gof_consistency( model_spline)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

ntrials = 100;
data = model_spline.data;
nelectrodes = size(data,1);
N = size(data,2);
trials_true     = zeros(nelectrodes,N,ntrials);
trials_spline     = zeros(nelectrodes,N,ntrials);


for i = 1:ntrials
    trials_spline(:,:,i) = simulate_data(model_spline);
    trials_true(:,:,i) = data;
end
s = ntrials*size(data,2);
Y=reshape(trials_spline,[nelectrodes s]);
X=reshape(trials_true,[nelectrodes s]);

Rr = (X*X')/(s-1);           % covariance estimate
Rs = (Y*Y')/(s-1);           % covariance estimate

cons = 1 - norm(Rs-Rr)/norm(Rr); % compare matrix norms


end

