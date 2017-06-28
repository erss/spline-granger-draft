
function b = three_node_sim_3
a1 = 0.07*[hann(20)', -0.5*ones(20,1)']';
    a2 = 0.03*[-0.5*ones(20,1)', hann(20)']';
    b = zeros(3,3,40);
    b(1,3,:) = a1;                             % Model coefficients
    b(2,3,:) = a1;
    b(3,3,:) = a2; 
end