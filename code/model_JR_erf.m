% MODEL_JR This function implements the state space representation of
% the nonlinear neural mass model equations.
% 
% Inputs: x - the current state
%         params - a structure of parameters
%         mode - 'transition' returns the new state,
%                'jacobian' returns the Jacobian at the current state.
%
% Jan 2013
%
function out = model_JR_erf(x,mode,params)

dt = params.dt;         % time step
mu = params.mu;        % this is just the mean when this function is used in estimation

% Get current parameters
%
e_0 = params.e_0;        % maximal firing parameter
v_0 = params.v_0;        % firing threshold
r = params.r;              % sigmoid slope

A = params.A;          % excitatory gain
a = params.a;          % excitatory time constant
B = params.B;         % inhibitory gain
b = params.b;         % inhibitory time constant

C1 = params.C1;        % connectivity parameter - excitatory feedback PSP
C2 = params.C2;        % connectivity parameter - excitatory feedback firing
C3 = params.C3;        % connectivity parameter - inhibitory feedback PSP
C4 = params.C4;        % connectivity parameter - excitatory feedback firing


% Linear component of model
%
F = [1, dt, 0, 0, 0, 0; ...
    -a^2*dt, 1-2*a*dt, 0, 0, 0, 0; ...
    0, 0, 1, dt, 0, 0; ...
    0, 0, -a^2*dt, 1-2*a*dt, 0, 0; ...
    0, 0, 0, 0, 1, dt; ...
    0, 0, 0, 0, -b^2*dt, 1-2*b*dt];

% Calculate the width of equivalent sigmoid using error function. See
% test_sigmoid_erf.m function for optimisation results
%
rErf = 1.699/r;  

% Sigmoid functions using error function
%
f_p = 0.5 * (1 + erf( ((x(3)-x(5)) - v_0)/(sqrt(2).*rErf) )) ;   % pyramidal firing
f_e = 0.5 * (1 + erf( (C1*x(1) - v_0)/(sqrt(2).*rErf) ));       % excitatory feedback firing
f_i = 0.5 * (1 + erf( (C3*x(1) - v_0)/(sqrt(2).*rErf) ));       % inhibitory feedback firing


if mode(1)=='t';
    
    % Nonlinear component
    %
    gx = [0; ...
        dt*A*a*2*e_0*f_p; ...
        0; ...
        dt*A*a*C2*2*e_0*f_e; ...
        0 ; ...
        dt*B*b*C4*2*e_0*f_i];
    
    % Constant component
    %
    c = [0; 0; 0; dt*A*a*mu; 0; 0];
    
    % Nonlinear transition model
    %
    out = F*x + gx + c;
else
    % Linearise g()
    %
    f_p_derivative = 2*e_0*1/sqrt(2*pi)*exp(-(((x(3)-x(5)) - v_0)/(sqrt(2).*rErf)).^2)./rErf;         % pyramidal firing
    f_e_derivative = 2*e_0*C1/sqrt(2*pi)*exp(-((C1*x(1) - v_0)/(sqrt(2).*rErf)).^2)./rErf;      % excitatory feedback firing
    f_i_derivative = 2*e_0*C3/sqrt(2*pi)*exp(-((C3*x(1) - v_0)/(sqrt(2).*rErf)).^2)./rErf;      % inhibitory feedback firing

    G = [0, 0, 0, 0, 0, 0; ...
        0,0,dt*A*a*f_p_derivative,0,-dt*A*a*f_p_derivative, 0; ...
        0, 0, 0, 0, 0, 0; ...
        dt*A*a*C2*f_e_derivative, 0, 0, 0, 0, 0; ...
        0, 0, 0, 0, 0, 0; ...
        dt*B*b*C4*f_i_derivative, 0, 0, 0, 0, 0];
    
    % Jacobian
    %
    out = F + G;
end


% get current states
% v_f = x(1);
% z_f = x(2);
% v_p1 = x(3);
% z_p1 = x(4);
% v_p2 = x(5);
% z_p2 = x(6);

% state update equation - propagate states
% x(1) = v_f + dt*z_f;
% x(2) = z_f + dt*(A*a*f_p - 2*a*z_f - a^2*v_f);
% x(3) = v_p1 + dt*z_p1;
% x(4) = z_p1 + dt*(A*a*(mu + C2*f_e) - 2*a*z_p1 - a^2*v_p1);
% x(5) = v_p2 + dt*z_p2;
% x(6) = z_p2 + dt*(B*b*C4*f_i - 2*b*z_p2 - b^2*v_p2);