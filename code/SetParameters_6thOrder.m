

% set parameters for the 6th order JR model

% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%% this set is for alpha rhythm
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
mu = 220;        % this is just the mean when this function is used in estimation
sigma = 5.74;


% get current parameters
e_0 = 2.5;        % maximal firing parameter
v_0 = 6;        % firing threshold
r = 0.56;              % sigmoid slope

A = 3.25;          % excitatory gain
a = 100;          % excitatory time constant
B = 22;         % inhibitory gain
b = 50;         % inhibitory time constant

% connectivity parameters from J and R model
C = 135; % 270;        % 135 for alpha activity
C1 = C;        % connectivity parameter - excitatory feedback PSP
C2 = 0.8*C;        % connectivity parameter - excitatory feedback firing
C3 = 0.25*C;        % connectivity parameter - inhibitory feedback PSP
C4 = 0.25*C;        % connectivity parameter - excitatory feedback firing


% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% %% this set is for seizure
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
mu = 90;        % this is just the mean when this function is used in estimation
sigma = 1.74;


% get current parameters
e_0 = 2.5;        % maximal firing parameter
v_0 = 6;        % firing threshold
r = 0.56;              % sigmoid slope

A = 8.25;          % excitatory gain
a = 100;          % excitatory time constant
B = 22;         % inhibitory gain
b = 50;         % inhibitory time constant

% connectivity parameters from J and R model
C = 135; % 270;        % 135 for alpha activity
C1 = C;        % connectivity parameter - excitatory feedback PSP
C2 = 0.8*C;        % connectivity parameter - excitatory feedback firing
C3 = 0.25*C;        % connectivity parameter - inhibitory feedback PSP
C4 = 0.25*C;        % connectivity parameter - excitatory feedback firing