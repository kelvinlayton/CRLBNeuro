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
delta = 0.1;  % Sampling period
t=0:delta:4;
N=length(t);

% Transition model (linear)
%
Fmat = [1.5, delta; 0 1.5];
f = @(x)(Fmat*x);
F = @(x)Fmat;

% Observation matrix (only observe position)
%
Hmat = [1 0];
H = @(x)Hmat;
% Define noise vectors
%
sigma_v = 1;  % Process noise
sigma_w = 1;   % Measurement noise

Q = sigma_v^2*[1 0; 0 1];   % Non-singular Q - either recursion works
% Q = sigma_v^2*[1 0; 0 0];   % Singular Q - must use Bergman recursion

R = sigma_w^2 .* eye(size(Hmat,1));

% Prior distribution
%
x0 = [3; 0];
m0 = x0;
P0 = diag(100*[(pi)^2 (pi)^2]);


%% 
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Compute the posterior Cramer-Rao bound (PCRB)
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
M = 1000;

crb1 = compute_pcrb_J(t,f,F,H,Q,R,m0,P0,M);
crb2 = compute_pcrb_P(t,f,F,H,Q,R,m0,P0,M);

norm(crb1-crb2)

figure
semilogy(t,crb1,'-o');
hold on;
semilogy(t,crb2,'-x');


return

%% 
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Test algorithm details
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


% Initialise variables
%
pcrb=zeros(2,N);
pcrb2=zeros(2,N);
J=zeros(2,2,N);
Pprior=zeros(2,2,N);
Pprior(:,:,1)=P0;

Pprior2=zeros(2,2,N);
Pprior2(:,:,1)=P0;
J(:,:,1)=inv(P0);

% Compute the PCRB using closed form expressions
%
U = Fmat'*inv(Q)*Fmat;
V = -Fmat'*inv(Q);
W = inv(Q) + Hmat'*inv(R)*Hmat;

Rinv = Hmat'*inv(R)*Hmat;
G = eye(2);


Ppost=zeros(2,2,N);
Ppost(:,:,1)=P0;  % P_{0|0}

for k=2:N
    
    % Recursively compute the Fisher information matrix
    %
    J(:,:,k) = W - V'*inv(J(:,:,k-1) + U)*V;
   
    % P_{1|0}
    Pprior(:,:,k) = Fmat*Ppost(:,:,k-1)*Fmat' + Q;
    
    % P_{1|1}
    Ppost(:,:,k) = inv( inv(Pprior(:,:,k)) + Rinv);
    
    % Compute the PCRB at the curret time
    %
    pcrb(:,k) = diag(inv(J(:,:,k)));
    pcrb2(:,k) = diag(Ppost(:,:,k));

end
% pcrb2(:,end)=[];
norm(pcrb-pcrb2)./norm(pcrb2)

semilogy(t,pcrb);
hold on;
semilogy(t,pcrb2,'x');
semilogy(t,abs(pcrb-pcrb2),'r');
% semilogy(t,abs(pcrb-pcrb3),'.-');