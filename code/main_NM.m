% This code repricates the structure of the Voss / Schiff FN simulation /
% estimation

%%
% clc
% close all
% clear


params = SetParametersNM('alpha');

params.dt = 0.001;


N = 5000;             	% number of samples
dT = params.dt;          	% sampling time step (global)
dt = 1*dT;            	% integration time step
nn = fix(dT/dt);      	% (used in for loop for forward modelling) the integration time step can be small that the sampling (fix round towards zero)
t = 0:dt:(N-1)*dt;



% Transition model
NStates = 4;                           
f = @(x)model_NM(x,'transition',params);
F = @(x)model_NM(x,'jacobian',params);

% Initialise trajectory state
x0 = zeros(NStates,1);                   % initial state
x = zeros(NStates,N); 
x(:,1) = mvnrnd(x0,10^2*eye(NStates));


% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Generate trajectory
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

% Euler integration
%
for n=1:N-1
    x(:,n+1) = f(x(:,n));
end

% Calculate noise covariance based on trajectory variance over time??
%   Why noise on all states?
%
Q = 10^2.*diag((0.4*std(x,[],2)*sqrt(dt)).^2);

% Initialise random number generator for repeatability
rng(0);

v = mvnrnd(zeros(NStates,1),Q,N)';

% Euler-Maruyama integration
for n=1:N-1
    x(:,n+1) = f(x(:,n)) + v(:,n);
end

H = [1 0 0 0];           % observation function
R = 1^2*eye(1);

w = mvnrnd(zeros(size(H,1),1),R,N)';
y = H*x + w;

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
plot(t,x([1],:)'); hold on;
plot(t,m([1],:)','--');
ax2=subplot(212);
plot(t,y)

linkaxes([ax1 ax2],'x');

