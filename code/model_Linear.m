% MODEL_Linear This function implements the state space representation of
% a simple linear equation
% 
% Inputs: x - the current state
%         params - a structure of parameters
%         mode - 'transition' returns the new state,
%                'jacobian' returns the Jacobian at the current state.
%
% Jan 2013
%
function out = model_Linear(x,mode,params)

F = [1, params.dt 0; 0 1 0; 0 params.dt 1];
% F = [1, params.dt ; 0 1];

if mode(1)=='t';
    out = F*x;
else
    out = F;
end