% COMPUTE_PCRB_J Compute the Posterior CRB using Tichavsky iteration
% 
% Inputs: t - vector of time points
%         f - transition function. For a regular Kalman fiter use 
%                 @(x)(F*x), where F is the transition matrix
%         F - transition matrix function (a function that takes the
%                 current state and returns the Jacobian). For a regular
%                 Kalman fiter use @(x)F, where F is the transition matrix
%         H - observation matrix function. For linear measurement 
%                 use @(x)H, where H is observation matrix
%         Q - process covariance
%         R - measurement covariance
%         m0 - mean of prior distribution
%         P0 - covariance of prior distribution
%         M - number of Monte Carlo samples
%
% Outputs: pcrb - the posterior CRB
%
% Kelvin Layton
% Feb 2013
%
function pcrb = compute_pcrb_J(t,f,F,H,Q,R,m0,P0,M)

N = length(t);
Nstates=length(m0);

% Initialise variables
%

J=zeros(Nstates,Nstates,N);
J(:,:,1)=inv(P0);
pcrb=zeros(Nstates,N);

% Initalise all trajectories
%
xk=mvnrnd(m0,P0,M)';


Rinv=inv(R);
Qinv=inv(Q);


% Compute the PCRB using a Monte Carlo approximation
%
for k=2:N
    U = zeros(Nstates);
    V = zeros(Nstates);
    W = zeros(Nstates);
    
    v = mvnrnd(zeros(Nstates,1),Q,M)';
    parfor i=1:M
        % Linearise transition function
        %
        Fhat = F(xk(:,i));
        
        % Sample the next time point for the current trajectory realisation
        %
        xk(:,i) = f(xk(:,i)) + v(:,i);

        % Linearise measurment function
        %
        Hhat = H(xk(:,i));
        
        U = U + Fhat'*Qinv*Fhat;
        V = V - Fhat'*Qinv;
        W = W + Qinv + Hhat'*Rinv*Hhat;

    end
    U=U./M;
    V=V./M;
    W=W./M;

    % Recursively compute the Fisher information matrix
    %
    J(:,:,k) = W - V'*inv(J(:,:,k-1) + U)*V;
    
    % Compute the PCRB at the current time
    %
    pcrb(:,k) = diag(inv(J(:,:,k)));

end