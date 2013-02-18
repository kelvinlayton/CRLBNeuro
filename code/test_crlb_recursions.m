% This script performs a quick validation of the Bergman recursion as
% equivalent to the Tichavsky recursion for non-singular Q matrices
%
% The advantage of the Bergman recursion is apparent for singular Q
% matrices
%
%


%%
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Simulate linear motion model and KF estimates
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

% close all

% Define time
%
delta = 0.01;  % Sampling period
N=30;
t=(0:N-1)*delta

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
sigma_v = 1;  % Process noise
sigma_w = 1;   % Measurement noise

Q = sigma_v^2*[1 0; 0 1];   % Non-singular Q - either recursion works
% Q = sigma_v^2*[1 0; 0 0];   % Singular Q - must use Bergman recursion

R = sigma_w^2 .* eye(size(H,1));

% Prior distribution
%
x0 = [3; 3];
m0 = x0;
P0 = diag(100*[(pi)^2 (pi)^2]);


%% 
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Compute the posterior Cramer-Rao bound (PCRB)
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


% Initialise variables
%
pcrb=zeros(2,N);
pcrb2=zeros(2,N);
J=zeros(2,2,N);
J2=zeros(2,2,N);
P=zeros(2,2,N);
P2=zeros(2,2,N);
P3=zeros(2,2,N);
Pkf=zeros(2,2,N);
P(:,:,1)=P0;
P2(:,:,1)=P0;
P3(:,:,1)=P0;
Pkf(:,:,1)=P0;
J(:,:,1)=inv(P0);
J2(:,:,1)=inv(P0);

pcrb(:,1) = diag(inv(J(:,:,1)));
pcrb2(:,1) = diag(inv(J2(:,:,1)));
pcrb3(:,1) = diag(P(:,:,1));
pcrb4(:,1) = diag(Pkf(:,:,1));
pcrb5(:,1) = diag(P2(:,:,1));
pcrb6(:,1) = diag(P3(:,:,1));

% Compute the PCRB using closed form expressions
%
U = F'*inv(Q)*F;
V = -F'*inv(Q);
W = inv(Q) + H'*inv(R)*H;

Rinv = H'*inv(R)*H;
G = eye(2);

B=sqrtm(Q);

for k=2:30
    
    % Recursively compute the Fisher information matrix
    %
    J(:,:,k) = W - V'*inv(J(:,:,k-1) + U)*V;
    J2(:,:,k) = H'*inv(R)*H + inv(Q) - inv(Q)*F*inv(J2(:,:,k-1) + F'*inv(Q)*F)*F'*inv(Q);
    
    P(:,:,k) = F*inv(inv(P(:,:,k-1)) + Rinv)*F' + G*Q*inv(G);
    
    Ppredict = F*Pkf(:,:,k-1)*F' + Q;
    S = H*Ppredict*H' + R;
    K = Ppredict*H'*inv(S);
    Pkf(:,:,k) = (eye(2)-K*H)*Ppredict;
    
    % Bergman thesis
    Qtil=inv(Q);
    Stil=F*inv(Q);
    Vtil=F*inv(Q)*F';
    Rtil=H*inv(R)*H';
    P3(:,:,k) = Qtil + Rtil - Stil'*inv(inv(P(:,:,k-1)) + Vtil)*Stil;
%     
    % Bergman paper
    Qtil_inv = B'*H'*inv(R)*H*B + eye(2);
    Rtil_inv = F'*H'*inv(R)*H*F;
    Stil = (B'*H'*inv(R)*H*F)';
    Up = inv(inv(P2(:,:,k-1)) + Rtil_inv);
    Up_inv = inv(P2(:,:,k-1)) + Rtil_inv;
    Delta=Qtil_inv - Stil'*Up_inv*Stil;
    P2(:,:,k) = F*Up*F' + (F*Up*Stil-B)*inv(Delta)*(Stil'*Up*F'-B);
    
    % Compute the PCRB at the curret time
    %
    pcrb(:,k) = diag(inv(J(:,:,k)));
    pcrb2(:,k) = diag(inv(J2(:,:,k)));
    pcrb3(:,k) = diag(P(:,:,k));
    pcrb4(:,k) = diag(Pkf(:,:,k));
    pcrb5(:,k) = diag(P2(:,:,k));
    pcrb6(:,k) = diag(P3(:,:,k));
end

% plot(pcrb(:,1:30)','r');
hold on;
% plot(pcrb2(:,1:30)','rx');
plot(pcrb3(:,1:30)','o');
plot(pcrb4(:,1:30)','+');
% plot(pcrb5(:,1:30)','rs');
% plot(pcrb6(:,1:30)','rd');


return


%%
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Simulate linear motion model and KF estimates
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%% Compute the MSE of the extended Kalman filter 
%
num_trials = 10000;
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
    m = extended_kalman_filter(z,f,Ffunc,H,Q,R,m0,P0);

    % Accumulate the estimation error
    %
    error = error + (x-m).^2;
end

% Calculate the mean squared error
%
mse = error ./ num_trials;

%%
semilogy(pcrb4(:,1:N)');
hold on;
semilogy(pcrb3(:,1:N)','o');
semilogy(mse(:,1:N)','x')

