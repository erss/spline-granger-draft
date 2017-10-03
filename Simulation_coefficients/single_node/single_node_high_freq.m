
function b = single_node_high_freq
% Returns coefficients for univariate AR(2) model that generate a high
% frequency signal.
b = zeros(1,1,2);
b(1,1,:) = [0.9 -0.8];
end