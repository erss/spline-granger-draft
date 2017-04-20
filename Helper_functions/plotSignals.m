function plotSignals( data )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

nm = nanmean(data, 1);
data = bsxfun(@minus, data, nm);
plot(data(1,:))
hgt = 0;
hold on
for i = 2:size(data,1)
    hgt = max(data(i,:)) -min(data(i,:))+hgt;
    data(i,:) =  data(i,:)+ hgt.*ones(size(data(1,:)));
    plot(data(i,:))

    
end

end

