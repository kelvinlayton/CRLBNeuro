% This code repricates the structure of the Voss / Schiff FN simulation /
% estimation
clc
close all
clear

N = 500;             	% number of samples
dT = 0.001;          	% sampling time step (global)
dt = 1*dT;            	% integration time step
nn = fix(dT/dt);      	% (used in for loop for forward modelling) the integration time step can be small that the sampling (fix round towards zero)

t = 0:dt:(N-1)*dt;


% Initialise random number generator for repeatability
rng(0);

M = 100;    % Mumber of Monte Carlo trials for MSE calculation


NStatesModels = [6 6];  

% Model parameters
params{1} = SetParametersLopes('alpha');
params{2} = SetParametersLopes('alpha');


% Transition functions
models{1} = @model_Lopes;
models{2} = @model_Lopes;

% Measurement fnctions
Hmodels{1} = [1 0 -1 0 0 0];           % Lopes
Hmodels{2} = [1 0 -1 0 0 0];           % Lopes

% Indices to plot (only potentials, not derivatives)
displayStates{1} = [1 3 5];
displayStates{2} = [1 3 5];

noiseInd = [2 2 ];     % NEED TO CHECK THESE! (and adjust Wendling noise)
noiseScale = [1 1 ];

mseModels = cell(length(models),1);
crbModels = cell(length(models),1);

for iModel = 1:length(models)
    
parameters = params{iModel};
parameters.dt = dt;
A = parameters.A;
a = parameters.a;
sigma = parameters.sigma;

% Transition model
NStates = NStatesModels(iModel);                           
f = @(x)models{iModel}(x,'transition',parameters);
F = @(x)models{iModel}(x,'jacobian',parameters);

H = Hmodels{iModel};

displayInd = displayStates{iModel};

% Initialise trajectory state
x0 = zeros(NStates,1);                   % initial state

% Define input
Q = zeros(NStates);
Q(noiseInd(iModel),noiseInd(iModel)) = noiseScale(iModel)*(sqrt(dt)*A*a*sigma)^2;

% Prior distribution (defined by m0 & P0)
%
m0 = x0;
if iModel==1
%     P0 = 100.^2*eye(NStates);
    P0 = 50.^2*diag([1 0 1 0 1 0]) + 50^2*diag([0 1 0 1 0 1]);
else
%     P0 = 200.^2*eye(NStates);
    P0 = 5000.^2*diag([1 0 1 0 1 0]) + 5000^2*diag([0 1 0 1 0 1]);
end


Hfun = @(x)H;
R = 10^2*eye(1);
    
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Compute CRLB
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
crbAll = compute_pcrb_P(t,f,F,Hfun,Q,R,m0,P0,200);

crb = crbAll(:,end);


semilogy(t,crbAll(displayInd,:))
iModel
hold on;

crbModels{iModel} = crb(displayInd);

end

return;

save crb_convergence_results.mat crbModels

plotConvergenceResults
