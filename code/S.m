% Normalised sigmoid function
%
function out = S(v,v_0,r)

out = 1 ./ (1 + exp(r*(v_0 - v)));