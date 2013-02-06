
% this code repricates the structure of the Voss / Schiff FN simulation /
% estimation

clear all
close all
clc

tic

NStates = 6;                               % in the non-augmented JR model

N = 100000;                                % number of samples
dT = 0.001;                               % sampling time step (global)
dt = 1*dT;                                 % integration time step
nn = fix(dT/dt);                         % (used in for loop for forward modelling) the integration time step can be small that the sampling (fix round towards zero)

t = 0:dt:(N-1)*dt;

% preallocate for speed
x0 = zeros(NStates,N);              % true state trajectory

% Initial conditions
x0(:,1) = zeros(NStates,1);                   % initial state

% intial true parameters value
SetParameters_6thOrder;

parameters = {dt,...
    mu,...
    e_0,...
    v_0,...
    r,...
    A,...         % excitatory gain
    a,...          % excitatory time constant
    B,...        % inhibitory gain
    b,...        % inhibitory time constant
    C1,...       % connectivity parameter - excitatory feedback PSP
    C2,...       % connectivity parameter - excitatory feedback firing
    C3,...       % connectivity parameter - inhibitory feedback PSP
    C4};       % connectivity parameter - excitatory feedback firing

% define input
e = sqrt(dt)*A*a*sigma*randn(N,1);

% Euler-Maruyama integration
for n=1:N-1
    x0(:,n+1) = JRint_6thOrder(x0(:,n), parameters)  + [0; 0; 0; e(n); 0; 0];
end

C = [0 0 1 0 -1 0];           % observation function
y = C*x0;

figure
plot(t,x0')
plot(t,y);