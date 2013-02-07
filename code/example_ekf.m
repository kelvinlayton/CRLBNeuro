%% Simulate pendulum model and EKF estimates
%
close all

% Define time
%
delta = 0.01;  % Sampling period
t=0:delta:4;
N=length(t);

% Transition model
%
f = @(x)model_pendulum(x,delta,'transition');
F = @(x)model_pendulum(x,delta,'jacobian');

% Observation matrix (only observe angle)
%
H = [1 0];

% Define noise vectors
%
sigma_v = 0.01;  % Process noise
sigma_w = 0.1;   % Measurement noise
Q = sigma_v^2*[1 0; 0 1];
R = sigma_w^2 .* eye(size(H,1));

w = mvnrnd(zeros(size(R,1),1),R,N)';
v = mvnrnd(zeros(size(Q,1),1),Q,N)';

% Starting state
%
x0 = [pi/4; 0];

% Calculate state trajectory
%
x = zeros(length(x0),N);
x(:,1)=x0;
for i=2:N
    x(:,i) = f(x(:,i-1)) + v(:,i);
end

% Define observation vector
%
z = H*x + w;

% Prior distribution (defined by m0 & P0)
%
m0 = x0;
P0 = diag([(2*pi)^2 (2*pi)^2]);

% Apply EKF filter
%
m = extended_kalman_filter(z,f,F,H,Q,R,m0,P0);

% Plot trajecotry and EKF estimates
%
figure;
plot(t,x); hold on;
plot(t,m,'--');
xlabel('Time (s)');
ylabel('States');
grid on;
legend({'x_1 - angle','x_2 - velocity','x_1 estimate','x_2 estimate'})
% export_fig('../figures/pendulum_trajectory.pdf')


%% Compute the posterior Cramer-Rao bound (PCRB)
%

M = 100;    % Number of Monte Carlo samples

% Initialise variables
%

Qinv=inv(Q);
Rinv=inv(R);
J=zeros(2,2,N);
J(:,:,1)=inv(P0);
pcrb=zeros(2,N);

% Initalise all trajectories
%
xk=zeros(length(m0),M);
for i=1:M
    xk(:,i)=mvnrnd(m0,P0)';
end

% Compute the PCRB using a Monte Carlo approximation
%
for k=2:N
    U = zeros(2,2);
    V = zeros(2,2);
    
    v = mvnrnd([0 0]',Q,M)';
    parfor i=1:M
        % Sample the next time point for the current trajectory realisation
        %
        xk(:,i) = f(xk(:,i)) + v(:,i);

        % Compute the PCRB terms for the current trajectory realisation
        %
        Fhat = F(xk(:,i));
        U = U + Fhat'*Qinv*Fhat;
        V = V + -Fhat'*Qinv;
    end
    U=U./M;
    V=V./M;
    W=Qinv + H'*Rinv*H;
    
    % Recursively compute the Fisher information matrix
    %
    J(:,:,k) = W - V'*inv(J(:,:,k-1) + U)*V;

    % Compute the PCRB at the curret time
    %
    pcrb(:,k) = diag(inv(J(:,:,k)));
end


%% Compute the MSE of the extended Kalman filter 
%
num_trials = 1000;
error = zeros(2,N);

parfor r=1:num_trials
    
    % Create new trajectory realisation
    %
    v = mvnrnd([0 0]',Q,N)';
    x = zeros(2,N);
    x(:,1)=x0; %mvnrnd(m0,P0)';
    for i=2:N
        x(:,i) = f(x(:,i-1)) + v(:,i);
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
semilogy(t,pcrb,'-');
grid on;
legend({'MSE x_1','MSE x_2','PCRB x_1','PCRB x_2'})
xlabel('Time (s)');
ylabel('MSE');
% export_fig('../figures/pendulum_pcrb.pdf')