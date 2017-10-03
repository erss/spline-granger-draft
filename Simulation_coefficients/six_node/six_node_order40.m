
function b = six_node_order40
% Returns coefficients for multivariate autoregressive model with 6 signals:
% . b is 6x6x40 such that entry b(i,j,:) contain the model coefficients for
%     signal j's influence on signal i.
a1 = 0.07*[0.1*hann(20)', -0.5*ones(20,1)']';
a2 = 0.03*[-0.5*ones(20,1)', hann(20)']';    
%a3 = -0.3*ones(size(a1));
a3 = -.3*sin(1:40)./(1:40); 
%a3 = -.3*sin(.5:.5:20)./(.5:.5:20);
a4 = 0.02*[hann(30)', -1*ones(1,10)]';
%a4 = 0.02*[hann(30)', -0.2*ones(1,10)]';
a5 = 0.1*[-0.8*ones(1,30), -0.8*hann(10)']';
a6 = 0.25*[-0.6*ones(20,1)', -0.2*hann(20)']'; 

b = zeros(6,6,40);
b(1,1,:) = a1;
b(1,2,:) = a2;
b(2,2,:) = a1;
b(2,4,:) = a4;
b(3,3,:) = a3;
b(3,4,:) = a2;
b(4,4,:) = a4;
b(5,2,:) = a5;
b(5,5,:) = a4;
b(6,6,:) = a6;
b(6,3,:) = a4;
end