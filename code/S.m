function out = S(v)

e_0 = 2.5;
r = 0.56;
v_0 = 6;

out = 2*e_0 ./ (1 + exp(r*(v_0 - v)));