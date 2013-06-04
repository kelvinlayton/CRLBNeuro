
clc
close all
clear

N = 400;             	% number of samples
dT = 0.001;          	% sampling time step (global)
dt = 1*dT;            	% integration time step
nn = fix(dT/dt);      	% (used in for loop for forward modelling) the integration time step can be small that the sampling (fix round towards zero)

t = 0:dt:(N-1)*dt;

% Intial true parameter values

paramsL = SetParametersLopes('alpha');
paramsJR = SetParametersJR('seizure');
paramsW = SetParametersWendling('spikes');

paramsL.dt = dt;
paramsJR.dt = dt;
paramsW.dt = dt;

% Initialise random number generator for repeatability
rng(0);


A = paramsL.A;
a = paramsL.a;
sigma = paramsL.sigma;

Q_L = zeros(6);
Q_L(2,2) = (sqrt(dt)*A*a*sigma)^2;

Q_JR = zeros(6);
Q_JR(4,4) = (sqrt(dt)*A*a*sigma)^2;

Q_W = zeros(10);
Q_W(4,4) = (sqrt(dt)*A*a*sigma)^2;

H_L = [1 0 -1 0 0 0];
H_JR = [0 0 1 0 -1 0]; 
H_W = [0 0 1 0 -1 0 -1 0 0 0];           % observation function

R = 1^2*eye(1);


% Transition models
f_L = @(x)model_Lopes(x,'transition',paramsL);
F_L = @(x)model_Lopes(x,'jacobian',paramsL);

f_JR = @(x)model_JR(x,'transition',paramsJR);
F_JR = @(x)model_JR(x,'jacobian',paramsJR);

f_W = @(x)model_Wendling(x,'transition',paramsW);
F_W = @(x)model_Wendling(x,'jacobian',paramsW);

% Priors

m0_LJR = zeros(6,1);
P0_LJR = 100.^2*eye(6);

m0_W = zeros(10,1);
P0_W = 100.^2*eye(10);

%% 
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Compute the posterior Cramer-Rao bound (PCRB)
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

M = 500;    % Number of Monte Carlo samples

pcrb_L = compute_pcrb_P(t,f_L,F_L,@(x)H_L,Q_L,R,m0_LJR,P0_LJR,M);
pcrb_JR = compute_pcrb_P(t,f_JR,F_JR,@(x)H_JR,Q_JR,R,m0_LJR,P0_LJR,M);
pcrb_W = compute_pcrb_P(t,f_W,F_W,@(x)H_W,Q_W,R,m0_W,P0_W,M);

%%
l1=semilogy(t,sqrt(pcrb_L([1 3 5],:)),'-'); hold on;
l2=semilogy(t,sqrt(pcrb_JR([1 3 5],:)),'--');
l3=semilogy(t,sqrt(pcrb_W([1 3 5 7 9],:)),':');

