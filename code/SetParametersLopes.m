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
        params.mu = 220;        % this is just the mean when this function is used in estimation
        params.sigma = 5.74;
        
        % get current parameters
        params.e_0 = 2.5;        % maximal firing parameter
        params.v_0 = 6;        % firing threshold
        params.r = 0.56;              % sigmoid slope
        
        params.A = 3.25;          % excitatory gain
        params.a = 100;          % excitatory time constant
        params.B = 22;         % inhibitory gain
        params.b = 50;         % inhibitory time constant
        
        % connectivity parameters from Lopes model
        params.C1 = 100;        % connectivity parameter - excitatory feedback PSP
        params.C2 = 100;        % connectivity parameter - excitatory feedback firing

        
end
