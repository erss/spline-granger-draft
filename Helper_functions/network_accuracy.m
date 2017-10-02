function [ frac, tp, tn, fn, fp] = network_accuracy( A,B )
% NETWORK_ACCURACY computes the proportion of correctly identified edges
% and non edges. 
%
% INPUTS:
%  A,B= square binary matrices of same size
%
% OUTPUTS:
%  frac = fraction of correctly identified edges and non-edges
%   tp  = true positives
%   tn  = true negatives
%   fn  = false negatives
%   fp  = false postives
diff = A-B;
frac = length(find(diff~=0))/size(A,1)^2;
frac = 1-frac;

% Assume A is 'true' network & B is estimated network
% find edges
edges = find(A==1);
nonedges = find(A==0);

Be = B(edges);
tp = sum(Be == 1);
fn = sum(Be == 0);

Bne = B(nonedges);
tn = sum(Bne == 0);
fp = sum(Bne == 1);

end

