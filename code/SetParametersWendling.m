% Set parameters for the Wendling
%
function params = SetParametersWendling(mode)

if nargin<1
    mode = 'background';
end

% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% this set is for alpha rhythm
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


% set for 30 to 150 pps.
params.mu = 90;        % this is just the mean when this function is used in estimation
params.sigma = 5.74;%3.46;

% get current parameters
params.e_0 = 2.5;        % maximal firing parameter
params.v_0 = 6;        % firing threshold
params.r = 0.56;              % sigmoid slope

params.a = 100;          % excitatory gain
params.b = 50;         % inhibitory gain
params.g = 500;

% connectivity parameters from Wendling model
C = 135;
params.C1 = C;        % connectivity parameter - excitatory feedback PSP
params.C2 = 0.8*C;        % connectivity parameter - excitatory feedback firing
params.C3 = 0.25*C;
params.C4 = 0.25*C;
params.C5 = 0.3*C;
params.C6 = 0.1*C;
params.C7 = 0.8*C;
        
switch mode
    case 'background'
        params.A = 3.25;          % excitatory gain
        params.B = 22;         % inhibitory gain
        params.G = 10;

    case 'alpha'
        
        % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        % this set is for a seizure rhythm
        % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        params.A = 5;          % excitatory gain
        params.B = 15;         % inhibitory gain
        params.G = 4;
     
    case 'spikes'
        
        % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        % this set is for spikey bits before seizures
        % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        params.A = 5;          % excitatory gain
        params.B = 25;         % inhibitory gain
        params.G = 10;

 
end
