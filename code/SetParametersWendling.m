% Set parameters for the 6th order JR model
%
function params = SetParametersLopes(mode)

if nargin<1
    mode = 'alpha';
end

switch mode
    case 'alpha'
        
        % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        % this set is for alpha rhythm
        % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        
        
        % set for 30 to 60 pps.
        params.mu = 80;        % this is just the mean when this function is used in estimation
        params.sigma = 5.74;
        
        % get current parameters
        params.e_0 = 2.5;        % maximal firing parameter
        params.v_0 = 6;        % firing threshold
        params.r = 0.56;              % sigmoid slope
        
        params.A = 3.25;          % excitatory gain
        params.a = 100;          % excitatory time constant
        params.B = 22;         % inhibitory gain
        params.b = 50;         % inhibitory time constant
        params.G = 10;
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
        
end
