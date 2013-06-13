% This code repricates the structure of the Voss / Schiff FN simulation /
% estimation
clc
close all
clear

N = 2000;             	% number of samples
dT = 0.001;          	% sampling time step (global)
dt = 1*dT;            	% integration time step
nn = fix(dT/dt);      	% (used in for loop for forward modelling) the integration time step can be small that the sampling (fix round towards zero)

t = 0:dt:(N-1)*dt;

% Intial true parameter values
% mode = 'background';
mode = 'alpha';       % this is more like alpha!!!
% mode = 'spikes';
parameters = SetParametersWendling(mode);
parameters.dt = dt;
A = parameters.A;
a = parameters.a;
sigma = parameters.sigma;

% Initialise random number generator for repeatability
rng(0);

% Transition model
NStates = 10;                           
f = @(x)model_Wendling(x,'transition',parameters);
F = @(x)model_Wendling(x,'jacobian',parameters);

% Initialise trajectory state
x0 = zeros(NStates,1);                   % initial state
x = zeros(NStates,N); 
x(:,1) = x0;


% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Generate trajectory
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

% Define input
Q = zeros(NStates);
Q(4,4) = 0.4225*(sqrt(dt)*A*a*sigma)^2;         % note! the scaling makes the noise the same as the other models
v = mvnrnd(zeros(NStates,1),Q,N)';

% Euler-Maruyama integration

for n=1:N-1
    x(:,n+1) = f(x(:,n))  + v(:,n);
end

H = [0 0 1 0 -1 0 -1 0 0 0];           % observation function
R = 1^2*eye(1);
w = mvnrnd(zeros(size(H,1),1),R,N)';
y = H*x + w;

% save(['Wendling_' mode '.mat'],'y','t') 
figure
plot(t,-H*x)

%% Run EKF for this model

% Prior distribution (defined by m0 & P0)
%
m0 = x0;
P0 = 100.^2*eye(NStates);


% Apply EKF filter
%
m = extended_kalman_filter(y,f,F,H,Q,R,m0,P0);

figure
ax1=subplot(211);
plot(t,x([3 5 7],:)'); hold on;
plot(t,m([3 5 7],:)','--');
ax2=subplot(212);
plot(t,y)

linkaxes([ax1 ax2],'x');
