%% Simulate pendulum model and EKF estimates
%
close all

% Define time
%
delta = 0.001;  % Sampling period
N=400;
t=(0:N-1)*delta;

% Transition model
%
f = @(x)model_pendulum(x,delta,'transition');
F = @(x)model_pendulum(x,delta,'jacobian');

% Observation matrix (only observe angle)
%
H = [1 0];

% Define noise covariances
%
Q = (0.01)^2*[1 0; 0 1];         % Process noise
R = (0.1)^2 .* eye(size(H,1));    % Measurement noise

w = mvnrnd(zeros(size(R,1),1),R,N)';
v = mvnrnd(zeros(size(Q,1),1),Q,N)';

% Starting state
%
x0 = [pi/2; 0];

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
P0 = diag(0.001*[(2*pi)^2 (2*pi)^2]);

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

pcrb = compute_pcrb_P(t,f,F,@(x)H,Q,R,m0,P0,M);


%% Compute the MSE of the extended Kalman filter 
%
num_trials = 100;
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
semilogy(t,sum(mse),'x-')
hold on;
semilogy(t,sum(pcrb),'o-');
grid on;
legend({'MSE','PCRB'})
xlabel('Time (s)');
ylabel('MSE');
% export_fig('../figures/pendulum_pcrb.pdf')

% figure
% fnorm(1)=nan;
% plot(t,fnorm);
% grid on;
% xlabel('Time (s)');
% ylabel('Norm of F');
