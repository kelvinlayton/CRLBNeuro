% MODEL_PENDULUM This function implements the state space representation of
% the nonlinear pendulum equations.
% 
% Inputs: x - the current state
%         delta - the time step
%         mode - 'transition' returns the new state,
%                'jacobian' returns the Jacobian at the current state.
%
% Kelvin Layton
% Jan 2013
%
function [out] = model_pendulum(x, delta, mode)

    if mode(1)=='t';
        % Nonlinear transition model
        out = [x(1) + x(2)*delta; ...
                x(2) - 100*sin(x(1))*delta];

    else
        % Jacobian
        out = [1, delta; ...
                -100*cos(x(1))*delta, 1];
    end
end