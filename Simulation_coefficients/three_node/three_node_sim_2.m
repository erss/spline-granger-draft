
function b = three_node_sim_2
a1 = 0.07*[hann(20)', -0.5*ones(20,1)']';
    a2 = 0.03*[-0.5*ones(20,1)', hann(20)']';
    a3 = -.3*sin(1:40)./(1:40);   
    b = zeros(3,3,40);                         % Model coefficients
    b(1,1,:) = a1;
    b(1,2,:) = a2;
    b(2,1,:) = a1;
    b(2,2,:) = a2;
    b(3,3,:) = a3;
  
end