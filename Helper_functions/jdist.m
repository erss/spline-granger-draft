function [ dj] = jdist( A, B)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

M00 = length(find(A == 0 & B==0));
M01 = length(find(A == 0 & B==1));
M10 = length(find(A == 1 & B==0));
M11 = length(find(A == 1 & B==1));

%M = [M00 M10; M01 M11];

dj = (M01 + M10)/(M01+M10+M11);

% X = [A(:)';B(:)'];
% dj2 = pdist(X,'jaccard');

end

