%% Simulate linear motion model and KF estimates
%
close all

% Define time
%
delta = 0.01;  % Sampling period
t=0:delta:4;
N=length(t);

% Transition model (linear)
%
F = [1, delta; 0 1];
f = @(x)(F*x);
Ffunc = @(x)F;

% Observation matrix (only observe position)
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
x0 = [3; 0];

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
m = extended_kalman_filter(z,Ffunc,H,Q,R,m0,P0);

% Plot trajecotry and EKF estimates
%
figure;
plot(t,x); hold on;
plot(t,m,'--');
xlabel('Time (s)');
ylabel('States');
grid on;
legend({'x_1 - position','x_2 - velocity','x_1 estimate','x_2 estimate'})
% export_fig('../figures/motion_trajectory.pdf')


%% Compute the posterior Cramer-Rao bound (PCRB)
%

% Initialise variables
%
pcrb=zeros(2,N);
J=zeros(2,2,N);
J(:,:,1)=inv(P0);

% Compute the PCRB using closed form expressions
%
U = F'*inv(Q)*F;
V = -F'*inv(Q);
W = inv(Q) + H'*inv(R)*H;

for k=2:N
    
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
    x(:,1)=mvnrnd(m0,P0)';
    for i=2:N
        x(:,i) = f(x(:,i-1)) + v(:,i);
    end

    % Generate new observations 
    %
    w = mvnrnd(zeros(1,size(H,1)),R,N)';
    z = H*x + w;

    % Apply EKF filter (really the KF in this case)
    %
    m = extended_kalman_filter(z,Ffunc,H,Q,R,m0,P0);

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
% export_fig('../figures/motion_pcrb.pdf')