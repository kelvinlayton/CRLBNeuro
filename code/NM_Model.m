% the input is a deterministic sinusoidal firing rate + a stochastic input

function x_tplus1 = NM_Model(x,params)

% the parameters
dt = params.dt;

e0 = params.e0;
r = params.r;
v0 = params.v0;

a = params.a;
b = params.b;

A = params.A;
B = params.B;

mu = params.mu;        % mean input firing rate.

% the states
v_e = x(1);
z_e = x(2);

v_i = x(3);
z_i = x(4);

% v_x = x(5)
% z_x = x(6);

C1 = 100;
C2 = 100;

% update firing rates
f_i = 2*e0 ./ (1 + exp(r*(v0 - (mu - v_i))));          % inhibitory population firing rate
f_e = 2*e0 ./ (1 + exp(r*(v0 - v_e)));    % excitatory population firing rate

% transition

% input from inbib
v_e_tplus1 = z_e*dt + v_e;
z_e_tplus1 = (B*b*C2*f_i - 2*b*z_e - b^2*v_e)*dt + z_e;

% input from excit
v_i_tplus1 = z_i*dt + v_i;
z_i_tplus1 = (A*a*C1*f_e - 2*a*z_i - a^2*v_i)*dt + z_i;

% external input
% v_x_tplus1 = z_x*dt + v_x;
% z_x_tplus1 = (A*a*mu - 2*a*z_x - a^2*v_x)*dt + z_x;


% % output
x_tplus1 = [v_e_tplus1 z_e_tplus1 v_i_tplus1 z_i_tplus1]';

% output with state for external input
% x_tplus1 = [v_e_tplus1 z_e_tplus1 v_i_tplus1 z_i_tplus1 v_x_tplus1 z_x_tplus1]';



