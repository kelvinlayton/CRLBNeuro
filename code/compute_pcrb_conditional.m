% COMPUTE_PCRB_P Compute the Posterior CRB using Bergman iteration
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
function pcrb = compute_pcrb_conditional(t,f,F,H,Q,R,m0,P0,theta,M)

N = length(t);
Nstates=length(m0);
Nparams=length(theta);
Nunknown=Nstates+Nparams;

% Initialise variables
%
P0all = blkdiag(P0,inf*ones(Nparams,Nparams));
P=zeros(Nunknown,Nunknown,N);
P(:,:,1)=P0all;

pcrb=zeros(Nunknown,N);


% Initalise all trajectories
%
xk=mvnrnd(m0,P0,M)';
xk=[xk; bsxfun(@times,theta,ones(1,M))];

Rinv=inv(R);



% Compute the PCRB using a Monte Carlo approximation
%
for k=2:N
    Fhat = zeros(Nunknown);
    Rinvhat = zeros(Nunknown);
    
    v = mvnrnd(zeros(Nunknown,1),Q,M)';
    parfor i=1:M
        % Sample the next time point for the current trajectory realisation
        %
        xk(:,i) = f(xk(:,i)) + v(:,i);

        % Compute the PCRB terms for the current trajectory realisation
        %
        Fhat = Fhat + F(xk(:,i));
        
        Hmat = H(xk(:,i));
        Rinvhat = Rinvhat + Hmat'*Rinv*Hmat;

    end
    Fhat=Fhat./M;
    Rinvhat=Rinvhat./M;
        
    % Recursively compute the Fisher information matrix
    %
    P(:,:,k) = inv( inv(Fhat*P(:,:,k-1)*Fhat' + Q) + Rinvhat);
    
    % Compute the PCRB at the current time
    %
    pcrb(:,k) = diag(P(:,:,k));
end