M=1000000;

x0=2;
P0=4;

v0=1;
Q0=1;

% Generate samples and propogate them through a Gaussian function. 
% (Gaussian is the derivative of the error function sigmoid)
%
x = mvnrnd(x0,P0,M);

f = 1./sqrt(2*pi*Q0).* exp(-(x-v0).^2./(2*Q0));

mean_sample = mean(f)

mean_theoretical = 1./sqrt(2*pi*(Q0+P0)) .* exp(-(x0-v0).^2./(2*(Q0+P0)))
