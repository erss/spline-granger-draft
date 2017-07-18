function [ mn sem] = confidencebds( data )
%UNTITLED9 Summary of this function goes here
%   Detailed explanation goes here
    mn  = mean(data);
    sem = 2*std(data)/sqrt(length(data));

end

