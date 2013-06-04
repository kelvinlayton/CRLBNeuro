% Derivative of normalised sigmoid function
%
function out = S_derivative(v,v_0,r)

out = r*S(v,v_0,r)*(1-S(v,v_0,r));