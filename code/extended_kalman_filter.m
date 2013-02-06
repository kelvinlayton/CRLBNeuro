% EXTENDED_KALMAN_FILTER Implements an extended Kalman filter (or regular 
% Kalman filter).
% 
% Inputs: z - measurements
%         Ffunc - transition matrix function (a function that takes the
%                 current state and returns the Jacobian). For a regular
%                 Kalman fiter use @(x)F, where F is the transition matrix
%         H - the observation matrix
%         Q - process covariance
%         R - measurement covariance
%         m0 - mean of prior distribution
%         P0 - covariance of prior distribution
%
% Outputs: m - the posterior mean (the estimated state)
%          P - the posterior covariance
%
% Kelvin Layton
% Jan 2013
%
function [m, P] = extended_kalman_filter(z,Ffunc,H,Q,R,m0,P0)

    % Initialise variables
    %
    N=length(z);
    m=zeros([length(m0),N]);
    P=zeros([size(P0),N]);

    % Run Kalman filter over the given data
    %
    m(:,1)=m0;
    P(:,:,1)=P0;
    for i=2:N
        
        % Get linearisation for EKF (note when Ffunc=@(x)F it returns F)
        %
        Fhat=Ffunc(m(:,i-1));
        
        % Prediction step
        %
        m(:,i) = Fhat*m(:,i-1);
        P(:,:,i) = Fhat*P(:,:,i-1)*Fhat' + Q;

        % Update step
        %
        K = P(:,:,i)*H'*inv((H*P(:,:,i)*H' + R));
        m(:,i) = m(:,i) + K*(z(:,i)-H*m(:,i));
        P(:,:,i) = (eye(length(m0))-K*H)*P(:,:,i);
    end
end