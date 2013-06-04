% the input is a deterministic sinusoidal firing rate + a stochastic input

function out = model_NM(x,mode,params)

% the parameters
dt = params.dt;

e_0 = params.e0;
r = params.r;
v0 = params.v0;

a = params.a;
b = params.b;

A = params.A;
B = params.B;

mu = params.mu;        % mean input firing rate.

C1 = 100;
C2 = 100;

% the states
v_e = x(1);
z_e = x(2);

v_i = x(3);
z_i = x(4);



% Linear component of model
%
F = [1, dt, 0, 0; ...
    -b^2*dt, 1-2*b*dt, 0, 0; ...
    0, 0, 1, dt; ...
    0, 0, -a^2*dt, 1-2*a*dt];

% Sigmoid functions
%
f_i = 1 ./ (1 + exp(r*(v0 - (mu - v_i))));     % inhibitory population firing rate
f_e = 1 ./ (1 + exp(r*(v0 - v_e)));            % excitatory population firing rate

if mode(1)=='t';
    
    % Nonlinear component
    %
    gx = [0; ...
        dt*B*b*C2*2*e_0*f_i; ...
        0; ...
        dt*A*a*C1*2*e_0*f_e];
    
    
    % Nonlinear transition model
    %
    out = F*x + gx ;
    
%     % input from inbib
%     v_e_tplus1 = z_e*dt + v_e;
%     z_e_tplus1 = (2*e_0*B*b*C2*f_i - 2*b*z_e - b^2*v_e)*dt + z_e;
% 
%     % input from excit
%     v_i_tplus1 = z_i*dt + v_i;
%     z_i_tplus1 = (2*e_0*A*a*C1*f_e - 2*a*z_i - a^2*v_i)*dt + z_i;
% 
% 
%     % % output
%     out2 = [v_e_tplus1 z_e_tplus1 v_i_tplus1 z_i_tplus1]';
%     assert(norm(out-out2)<1e-10)

else
    % Linearise g()
    %

    f_i_derivative = 2*e_0*r*f_i*(1-f_i);      % inhibitory feedback firing
    f_e_derivative = 2*e_0*r*f_e*(1-f_e);      % excitatory feedback firing

    G = [0, 0, 0, 0; ...
        0,0,-dt*B*b*C2*f_i_derivative,0; ...
        0, 0, 0, 0; ...
        dt*A*a*C1*f_e_derivative, 0, 0, 0];
    
    % Jacobian
    %
    out = F + G;
end





