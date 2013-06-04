% This code repricates the structure of the Voss / Schiff FN simulation /
% estimation

%%
clc
close all
clear

params.dt = 0.0001;

params.e0 = 2.5;
params.r = 0.56;
params.v0 = 6;

params.a = 100;
params.b = 50;

params.A = 3.25;
params.B = 22;


params.mu = 11;        % mean input mem potential.



N = 50000;             	% number of samples
dT = params.dt;          	% sampling time step (global)
dt = 1*dT;            	% integration time step
nn = fix(dT/dt);      	% (used in for loop for forward modelling) the integration time step can be small that the sampling (fix round towards zero)
t = 0:dt:(N-1)*dt;



% Transition model
NStates = 4;                           
f = @(x)NM_Model(x,params);

% Initialise trajectory state
x0 = zeros(NStates,1);                   % initial state
x = zeros(NStates,N); 
x(:,1) = x0;





% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Generate trajectory
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

% Euler integration
%
for n=1:N-1
    x(:,n+1) = f(x(:,n));
end

H = [1 0 0 0];           % observation function

figure
plot(H*x)

Q = zeros(NStates);
for n=1:NStates
    Q(n,n) = (0.4*std(x(n,:))*sqrt(dt))^2;
end
% Initialise random number generator for repeatability
rng(0);
v = mvnrnd(zeros(NStates,1),Q,N)';

% Euler-Maruyama integration
for n=1:N-1
    x(:,n+1) = f(x(:,n)) + v(:,n);
end

R = 1^2*eye(1);
w = mvnrnd(zeros(size(H,1),1),R,N)';
y = H*x + w;

figure
plot(y)
drawnow


%%

% initalise state estimate
x_hat = size(x);
x_hat(:,1) = mean(x,2);

