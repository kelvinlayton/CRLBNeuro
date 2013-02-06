function x = JRint_6thOrder(x,parameters)

dt = parameters{1};         % time step
mu = parameters{2};        % this is just the mean when this function is used in estimation

% get current parameters
e_0 = parameters{3};        % maximal firing parameter
v_0 = parameters{4};        % firing threshold
r = parameters{5};              % sigmoid slope

A = parameters{6};          % excitatory gain
a = parameters{7};          % excitatory time constant
B = parameters{8};         % inhibitory gain
b = parameters{9};         % inhibitory time constant

C1 = parameters{10};        % connectivity parameter - excitatory feedback PSP
C2 = parameters{11};        % connectivity parameter - excitatory feedback firing
C3 = parameters{12};        % connectivity parameter - inhibitory feedback PSP
C4 = parameters{13};        % connectivity parameter - excitatory feedback firing

% get current states
v_f = x(1);
z_f = x(2);
v_p1 = x(3);
z_p1 = x(4);
v_p2 = x(5);
z_p2 = x(6);

v_p = v_p1 - v_p2;                   % the minus is where the inhibition is introduced

f_p = 2*e_0 / (1 + exp(r*(v_0 - v_p)));           % pyramidal firing
f_e = 2*e_0 / (1 + exp(r*(v_0 - C1*v_f)));       % excitatory feedback firing
f_i = 2*e_0 / (1 + exp(r*(v_0 - C3*v_f)));       % inhibitory feedback firing

% state update equation - propagate states
x(1) = v_f + dt*z_f;
x(2) = z_f + dt*(A*a*f_p - 2*a*z_f - a^2*v_f);
x(3) = v_p1 + dt*z_p1;
x(4) = z_p1 + dt*(A*a*(mu + C2*f_e) - 2*a*z_p1 - a^2*v_p1);
x(5) = v_p2 + dt*z_p2;
x(6) = z_p2 + dt*(B*b*C4*f_i - 2*b*z_p2 - b^2*v_p2);
