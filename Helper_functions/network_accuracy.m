function [ frac ] = network_accuracy( A,B )
% FINDS PROPORTION OF CORRECTLY IDENTIFIED EDGES AND NON-EDGES

diff = A-B;
frac = length(find(diff~=0))/size(A,1)^2;
frac = 1-frac;

end

