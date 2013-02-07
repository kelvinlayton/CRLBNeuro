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
sigma_w = 1;   % Measurement noise

Q = sigma_v^2*[1 0; 0 1];   % Non-singular Q - either recursion works
% Q = sigma_v^2*[1 0; 0 0];   % Singular Q - must use Bergman recursion

R = sigma_w^2 .* eye(size(H,1));

% Prior distribution
%
x0 = [3; 0];
m0 = x0;
P0 = diag([(pi)^2 (pi)^2]);

%% 
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Compute the posterior Cramer-Rao bound (PCRB)
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


% Initialise variables
%
pcrb=zeros(2,N);
pcrb2=zeros(2,N);
J=zeros(2,2,N);
P=zeros(2,2,N);
P(:,:,1)=P0;
J(:,:,1)=inv(P0);

% Compute the PCRB using closed form expressions
%
U = F'*inv(Q)*F;
V = -F'*inv(Q);
W = inv(Q) + H'*inv(R)*H;

Rinv = H'*inv(R)*H;
G = eye(2);

for k=2:N
    
    % Recursively compute the Fisher information matrix
    %
    J(:,:,k) = W - V'*inv(J(:,:,k-1) + U)*V;
    P(:,:,k) = F*inv(inv(P(:,:,k-1)) + Rinv)*F' + G*Q*inv(G);
    % Compute the PCRB at the curret time
    %
    pcrb(:,k) = diag(inv(J(:,:,k)));
    pcrb2(:,k) = diag(P(:,:,k));
end

semilogy(t,pcrb);
hold on;
semilogy(t,pcrb2,'x');