function [ frac tp tn] = network_accuracy( A,B )
% FINDS PROPORTION OF CORRECTLY IDENTIFIED EDGES AND NON-EDGES

diff = A-B;
frac = length(find(diff~=0))/size(A,1)^2;
frac = 1-frac;

% Assume A is 'true' network & B is estimated network
% find edges
edges = find(A==1);
nonedges = find(A==0);

Be = B(edges);
tp = sum(Be == 1);

Bne = B(nonedges);
tn = sum(Bne == 0);

end

