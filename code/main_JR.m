% This code repricates the structure of the Voss / Schiff FN simulation /
% estimation
%


N = 1000;             	% number of samples
dT = 0.001;          	% sampling time step (global)
dt = 1*dT;            	% integration time step
nn = fix(dT/dt);      	% (used in for loop for forward modelling) the integration time step can be small that the sampling (fix round towards zero)

t = 0:dt:(N-1)*dt;

% Intial true parameter values
%
% mode='alpha';
mode='seizure';
parameters = SetParametersJR(mode);
parameters.dt = dt;
A=parameters.A;
a=parameters.a;
sigma=parameters.sigma;

% Initialise random number generator for repeatability
%
rng(0);

% Transition model
%
NStates = 6;                           
f = @(x)model_JR(x,'transition',parameters);
F = @(x)model_JR(x,'jacobian',parameters);

% Initialise trajectory state
%
x0 = zeros(NStates,1);                   % initial state
x = zeros(NStates,N); 
x(:,1) = x0;

QnonSingular = (sqrt(dt)*A*a*sigma)^2*eye(NStates);
Q = zeros(NStates);
Q(4,4) = (sqrt(dt)*A*a*sigma)^2;

R = 1^2*eye(1);

% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Generate trajectory
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

% Define input
%
v = mvnrnd(zeros(NStates,1),Q,N)';

% Euler-Maruyama integration
%
for n=1:N-1
    x(:,n+1) = f(x(:,n))  + v(:,n);
end

H = [0 0 1 0 -1 0];           % observation function
w = mvnrnd(zeros(size(H,1),1),R,N)';
y = H*x + w;

% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Run EKF
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

% Prior distribution (defined by m0 & P0)
%
m0 = x0;
P0 = F(x0)*QnonSingular*F(x0)';


% Apply EKF filter
%
m = extended_kalman_filter(y,f,F,H,Q,R,m0,P0);

figure
plot(t,x([1 3 5],:)'); hold on;
plot(t,m([1 3 5],:)','--');



%% 
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Compute the posterior Cramer-Rao bound (PCRB)
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

M = 1000;    % Number of Monte Carlo samples

% Initialise variables
%

P=zeros(NStates,NStates,N);
P(:,:,1)=P0;
pcrb=zeros(NStates,N);

% Initalise all trajectories
%
xk=zeros(NStates,M);
for i=1:M
    xk(:,i)=mvnrnd(m0,P0)';
end

% Closed form terms
%
Rinv = H'*inv(R)*H;


% Compute the PCRB using a Monte Carlo approximation
%
for k=2:N
    Fhat = zeros(NStates,NStates);
    
    v = mvnrnd(zeros(NStates,1),Q,M)';
    parfor i=1:M
        % Sample the next time point for the current trajectory realisation
        %
        xk(:,i) = f(xk(:,i)) + v(:,i);

        % Compute the PCRB terms for the current trajectory realisation
        %
        Fhat = Fhat + F(xk(:,i));
        
    end
    Fhat=Fhat./M;
        
    % Recursively compute the Fisher information matrix
    %
    P(:,:,k) = Fhat*inv(inv(P(:,:,k-1)) + Rinv)*Fhat' + Q;
    
    % Compute the PCRB at the current time
    %
    pcrb(:,k) = diag(P(:,:,k));
end



%% 
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Compute the MSE of the extended Kalman filter 
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

%
num_trials = 100;
error = zeros(NStates,N);

parfor r=1:num_trials
    
    % Create new trajectory realisation
    %
    v = mvnrnd(zeros(NStates,1),Q,N)';
    x = zeros(NStates,N);
    x(:,1)=mvnrnd(m0,P0)';
    for i=2:N
        x(:,i) = f(x(:,i-1)) +  v(:,i);
    end

    % Generate new observations 
    %
    w = mvnrnd(zeros(1,size(H,1)),R,N)';
    z = H*x + w;

    % Apply EKF filter
    %
    m = extended_kalman_filter(z,f,F,H,Q,R,m0,P0);

    % Accumulate the estimation error
    %
    error = error + (x-m).^2;
end

% Calculate the mean squared error
%
mse = error ./ num_trials;

%% Plot MSE and the PCRB
%
figure
semilogy(t,mse,'x-')
hold on;
semilogy(t,pcrb,'-','LineWidth',2);
title(mode);

grid on;
xlabel('Time (s)');
ylabel('MSE');
% export_fig(['../figures/nmm_pcrb_' mode '.pdf'])