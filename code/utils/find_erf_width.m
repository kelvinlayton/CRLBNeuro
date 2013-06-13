function rErf = find_erf_width(rSig)

v0Sig = 0;
v0Erf = 0;
x = linspace(-10,10,1000);

sig1 = 1 ./ (1 + exp(rSig*(v0Sig - x)));

rErf = lsqnonlin(@cost_func,5);


figure
plot(x,sig1)
hold
plot(x,0.5 .* ( 1 + erf((x - v0Erf)./(sqrt(2)*rErf))),'r')

function diff = cost_func(r)
    sig2 = 0.5 .* ( 1 + erf((x - v0Erf)./(sqrt(2)*r)) );
    diff = sig1 - sig2;
end

end