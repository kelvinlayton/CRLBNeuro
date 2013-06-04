% Set parameters for the 6th order JR model
%
function params = SetParametersNM(mode)

if nargin<1
    mode = 'alpha';
end

switch mode
    case 'alpha'
        
        params.e0 = 2.5;
        params.r = 0.56;
        params.v0 = 6;

        params.a = 100;
        params.b = 50;

        params.A = 3.25;
%         params.A = 5;
        params.B = 22;


        params.mu = 11;        % mean input mem potential.

        
end
