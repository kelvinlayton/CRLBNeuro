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

% states from input arguement
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

% Linear component of model
%
F = [1, dt, 0, 0, 0, 0, 0, 0, 0, 0; ...
     -a^2*dt, 1-2*a*dt, 0, 0, 0, 0, 0, 0, 0, 0; ...
     0, 0, 1, dt, 0, 0, 0, 0, 0, 0; ...
     0, 0, -a^2*dt, 1-2*a*dt, 0, 0, 0, 0, 0, 0; ...
     0, 0, 0, 0, 1, dt, 0, 0, 0, 0; ...
     0, 0, 0, 0, -b^2*dt, 1-2*b*dt, 0, 0, 0, 0; ...
     0, 0, 0, 0, 0, 0, 1, dt, 0, 0; ...
     0, 0, 0, 0, 0, 0, -g^2*dt, 1-2*g*dt, 0, 0; ...
     0, 0, 0, 0, 0, 0, 0, 0, 1, dt; ...
     0, 0, 0, 0, 0, 0, 0, 0, -b^2*dt, 1-2*b*dt; ...
     ];
 
Sig = @(v)S(v,v_0,r);
Sig_der = @(v)S_derivative(v,v_0,r);

% deterministic transition

if mode(1)=='t';
    
    % Nonlinear component
    %
    gx = [0; ...
          dt*A*a*2*e_0*Sig(y_1 - y_2 - y_3); ...
          0;
          dt*A*a*C2*2*e_0*Sig(C1*y_0); ...
          0; ...
          dt*B*b*C4*2*e_0*Sig(C3*y_0); ...
          0; ...
          dt*G*g*C7*2*e_0*Sig(C5*y_0 - C6*y_4); ...
          0; ...
          dt*B*b*2*e_0*Sig(C3*y_0)];
      
	% Constant component
    %
    c = [0, 0, 0, A*a*mu*dt, 0, 0, 0, 0, 0, 0]';
    
    % Nonlinear transition model
    %
    out = F*x + gx + c;
    
%     y_0_tplus1 = y_5*dt + y_0;              % EPSP induced from Py
%     y_5_tplus1 = (A*a*S(y_1 - y_2 - y_3) - 2*a*y_5 - a^2*y_0)*dt + y_5;
%     
%     y_1_tplus1 = y_6*dt + y_1;              % EPSP on Py from ex feedback
%     y_6_tplus1 = (A*a*(mu + C2*S(C1*y_0)) - 2*a*y_6 - a^2*y_1)*dt + y_6;
%     
%     y_2_tplus1 = y_7*dt + y_2;              % slow IPSP on Py
%     y_7_tplus1 = (B*b*C4*S(C3*y_0) - 2*b*y_7 - b^2*y_2)*dt + y_7;
%     
%     y_3_tplus1 = y_8*dt + y_3;              % fast IPSP on Py
%     y_8_tplus1 = (G*g*C7*S(C5*y_0 - C6*y_4) - 2*g*y_8 - g^2*y_3)*dt + y_8;
%     
%     y_4_tplus1 = y_9*dt + y_4;
%     y_9_tplus1 = (B*b*S(C3*y_0) - 2*b*y_9 - b^2*y_4)*dt + y_9;
%     
%     out2 = [y_0_tplus1 y_5_tplus1 y_1_tplus1 y_6_tplus1 y_2_tplus1 y_7_tplus1 y_3_tplus1 y_8_tplus1 y_4_tplus1 y_9_tplus1]';
%     assert(norm(out-out2)<1e-10)

else
    % Linearise g()
    %
    
    G = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0; ...
         0, 0, dt*A*a*2*e_0*Sig_der(y_1-y_2-y_3), 0, -dt*A*a*2*e_0*Sig_der(y_1-y_2-y_3), 0, -dt*A*a*2*e_0*Sig_der(y_1-y_2-y_3), 0, 0, 0; ...
         0, 0, 0, 0, 0, 0, 0, 0, 0, 0; ...
         dt*A*a*C2*2*e_0*C1*Sig_der(C1*y_0), 0, 0, 0, 0, 0, 0, 0, 0, 0; ...
         0, 0, 0, 0, 0, 0, 0, 0, 0, 0; ...
         dt*B*b*C4*2*e_0*C3*Sig_der(C3*y_0), 0, 0, 0, 0, 0, 0, 0, 0, 0; ...
         0, 0, 0, 0, 0, 0, 0, 0, 0, 0; ...
         dt*G*g*C7*2*e_0*C5*Sig_der(C5*y_0-C6*y_4), 0, 0, 0, 0, 0, 0, 0, -dt*G*g*C7*2*e_0*C6*Sig_der(C5*y_0-C6*y_4), 0; ...
         0, 0, 0, 0, 0, 0, 0, 0, 0, 0; ...
         dt*B*b*2*e_0*C3*Sig_der(C3*y_0), 0, 0, 0, 0, 0, 0, 0, 0, 0; ...
        ];
    
    % Jacobian
    %
    out = F + G;
end


