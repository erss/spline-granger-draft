
function b = single_node_low_freq
%  model_coefficients = 0.9;
model_coefficients = [0.3 0.3];
% model_coefficients = [0.9 -0.1]
         
b = zeros(1,1,length(model_coefficients));
b(1,1,:) = model_coefficients;

end