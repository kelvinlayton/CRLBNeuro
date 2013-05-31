%% Wendling simulation

function out = model_Wendling(x,mode,params)

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
G = params.G;
g = params.g;

C1 = params.C1;        % connectivity parameter - excitatory feedback PSP
C2 = params.C2;        % connectivity parameter - excitatory feedback firing
C3 = params.C3;        
C4 = params.C4;        
C5 = params.C5;        
C6 = params.C6;        
C7 = params.C7;

% states
y_0 = x(1);        % 
y_5 = x(2);        % derivative of the above

y_1 = x(3);        % 
y_6 = x(4);        % derivative of the above

y_2 = x(5);         % 
y_7 = x(6);         % derivative of above

y_3 = x(7);         % 
y_8 = x(8);         % derivative of above

y_4 = x(9);         % 
y_9 = x(10);         % derivative of above


f_v_p = 2*e_0 ./ (1 + exp(r*(v_0 - (y_1 - y_2 - y_3))));
f_v_i = 2*e_0 ./ (1 + exp(r*(v_0 - C3*y_0)));      % firing rate of excitatory population
f_v_e = 2*e_0 ./ (1 + exp(r*(v_0 - C1*y_0)));      % firing rate of inhibitory population
f_v_g = 2*e_0 ./ (1 + exp(r*(v_0 - (C5*y_0 - C6*y_4))));

if mode(1)=='t';
    
    y_0_tplus1 = y_5*dt + y_0;
    y_5_tplus1 = (A*a*f_v_p - 2*a*y_5 - a^2*y_0)*dt + y_5;
    
    y_1_tplus1 = y_6*dt + y_1;
    y_6_tplus1 = (A*a*(mu - C2*f_v_e) - 2*a*y_6 - a^2*y_1)*dt + y_6;
    
    y_2_tplus1 = y_7*dt + y_2;
    y_7_tplus1 = (B*b*C4*f_v_i - 2*b*y_7 - b^2*y_2)*dt + y_7;
    
    y_3_tplus1 = y_8*dt + y_3;
    y_8_tplus1 = (G*g*C7*f_v_g - 2*g*y_8 - g^2*y_3)*dt + y_8;
    
    y_4_tplus1 = y_9*dt + y_4;
    y_9_tplus1 = (B*b*f_v_i - 2*b*y_9 - b^2*y_4)*dt + y_9;
    
    out = [y_0_tplus1 y_5_tplus1 y_1_tplus1 y_6_tplus1 y_2_tplus1 y_7_tplus1 y_3_tplus1 y_8_tplus1 y_4_tplus1 y_9_tplus1]';
    
end

