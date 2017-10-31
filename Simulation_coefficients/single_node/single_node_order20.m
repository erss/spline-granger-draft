function b = single_node_order20
% Returns coefficients for univariate AR(20) model
b = zeros(1,1,20);
% b(1,1,:) = [-0.2314 %%%%% spectral peak, high freq
%     0.1613
%     0.1405
%     0.0741
%     -0.0098
%     -0.0836
%     -0.1193
%     -0.1048
%     -0.0585
%     0.0011
%     0.0559
%     0.0874
%     0.0837
%     0.0614
%     0.0289
%     -0.0050
%     -0.0316
%     -0.0423
%     -0.0285
%     0.0185];
% 
x = [1 2 5 10 15 20 21];
y = [-.02314 0.1613 -0.1193 0.0559 0.005 0.0185 0];
xx = 1:21;
yy = spline(x,y,xx);
plot(x,y,'o',xx,yy)
yy(5)= yy(6);
yy(6) = mean([yy(5) yy(7)]);
yy(7) =mean([yy(6) yy(8)]);
yy(4) = yy(6);
yy(3) = 0.05;
yy(2) = 0.1;
plot(1:21,yy,'o',xx,yy)

b(1,1,:) = yy(1:end-1);


end