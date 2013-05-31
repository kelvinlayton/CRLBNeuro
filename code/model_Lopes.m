%% Lopex da Silva / van Rotterdam simulation

function out = model_Lopes(x,mode,params)

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

% states
v_e1 = x(1);        % excitatory membrane voltage induced from the input
z_e1 = x(2);        % derivative of the above

v_e2 = x(3);        % excitatory membrane potential induced from inhibitory feedback
z_e2 = x(4);        % derivative of the above

v_i = x(5);         % inhibitory memrane voltage induced from excitatory input
z_i = x(6);         % derivative of above

v_e = v_e1 - v_e2;      % total membrane potential of excitatory population - just used in sigmoid

f_v_i = 2*e_0 ./ (1 + exp(r*(v_0 - v_i)));      % firing rate of excitatory population
f_v_e = 2*e_0 ./ (1 + exp(r*(v_0 - v_e)));      % firing rate of inhibitory population

if mode(1)=='t';
    
    v_e1_tplus1 = z_e1*dt + v_e1;
    z_e1_tplus1 = (-2*a*z_e1 - a^2*v_e1 + a*A*mu)*dt + z_e1;

    v_e2_tplus1 = z_e2*dt + v_e2;
    z_e2_tplus1 = (-2*b*z_e2 - b^2*v_e2 + b*B*C2*f_v_i)*dt + z_e2;

    v_i_tplus1 = z_i*dt + v_i;
    z_i_tplus1 = (-2*a*z_i - a^2*v_i + a*A*C1*f_v_e)*dt + z_i;

    out = [v_e1_tplus1 z_e1_tplus1 v_e2_tplus1 z_e2_tplus1 v_i_tplus1 z_i_tplus1]';
    
end


