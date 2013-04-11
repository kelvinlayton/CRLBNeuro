%% Simulate pendulum model and EKF estimates
%
close all

% Define time
%
delta = 0.001;  % Sampling period
N=40;
t=(0:N-1)*delta;

% Transition model
%
f = @(x)model_pendulum(x,delta,'transition');
F = @(x)model_pendulum(x,delta,'jacobian');
fa = @(x)model_pendulum_augmented(x,delta,'transition');
Fa = @(x)model_pendulum_augmented(x,delta,'jacobian');

% Observation matrix (only observe angle)
%
H = [1 0];
Ha = [1 0 0];

% Define noise covariances
%
Q = (0.01)^2*eye(2);
Qa = blkdiag(Q,0);         % Process noise
R = (0.1)^2 .* eye(size(Ha,1));    % Measurement noise 

w = mvnrnd(zeros(size(R,1),1),R,N)';
v = mvnrnd(zeros(size(Qa,1),1),Qa,N)';

% Deterministic parameter
%
g = 100;

% Starting state
%
x0 = [pi/4; 0; g];

% Calculate state trajectory
%
x = zeros(length(x0),N);
x(:,1)=x0;
for i=2:N
    x(:,i) = fa(x(:,i-1)) + v(:,i);
end

% Define observation vector
%
z = Ha*x + w;

% Prior distribution (defined by m0 & P0)
%
m0 = [pi/2; 0];
m0a = [m0; g];
P0 = diag([(0.01*pi)^2 (0.01*pi)^2]);
P0param = 1000^2;
P0a = blkdiag(P0,P0param);

% Apply EKF filter
%
m = extended_kalman_filter(z,fa,Fa,Ha,Qa,R,m0a,P0a);
m2 = extended_kalman_filter(z,f,F,H,Q,R,m0,P0);

% Plot trajecotry and EKF estimates
%
figure;
plot(t,x); hold on;
plot(t,m,'--');
plot(t,m2,'-.');
xlabel('Time (s)');
ylabel('States');
grid on;
legend({'x_1 - angle','x_2 - velocity','x_1 estimate','x_2 estimate'},'Location','Best')
% export_fig('../figures/pendulum_trajectory.pdf')


%% Compute the posterior Cramer-Rao bound (PCRB)
%

M = 10000;    % Number of Monte Carlo samples

pcrb = compute_pcrb_P(t,f,F,@(x)H,Q,R,m0,P0,M);
pcrb_a = compute_pcrb_P(t,fa,Fa,@(x)Ha,Qa,R,m0a,P0a,M);
% pcrb2 = compute_pcrb_conditional(t,fa,Fa,@(x)Ha,Qa,R,m0(1:2),P0(1:2,1:2),g,M);

return

%% Compute the MSE of the extended Kalman filter 
%
num_trials = 100;
error = zeros(length(x0),N);

parfor r=1:num_trials
    
    % Create new trajectory realisation
    %
    v = mvnrnd(zeros(size(Q,1),1),Q,N)';
    x = zeros(length(m0),N);
    x(:,1)=mvnrnd(m0,P0)';
    x(3,1)=g;   % Not random
    for i=2:N
        x(:,i) = fa(x(:,i-1)) + v(:,i);
    end

    % Generate new observations 
    %
    w = mvnrnd(zeros(1,size(Ha,1)),R,N)';
    z = Ha*x + w;

    % Apply EKF filter
    %
    m = extended_kalman_filter(z,fa,Fa,Ha,Q,R,m0,P0);

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
% semilogy(t,(mse),'x-')
semilogy(t,(pcrb),'o-'); hold on;
semilogy(t,(pcrb_a),'s-');
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
