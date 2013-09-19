% This code repricates the structure of the Voss / Schiff FN simulation /
% estimation
clc
close all
clear

N = 100;             	% number of samples
dT = 0.001;          	% sampling time step (global)
dt = 1*dT;            	% integration time step
nn = fix(dT/dt);      	% (used in for loop for forward modelling) the integration time step can be small that the sampling (fix round towards zero)

t = 0:dt:(N-1)*dt;


% Initialise random number generator for repeatability
rng(0);

M = 100;    % Mumber of Monte Carlo trials for MSE calculation


NStatesModels = [6 6 6];  

% Model parameters
params{1} = SetParametersLopes('alpha');
params{2} = SetParametersLopes('alpha');
params{3} = SetParametersLopes('alpha');


% Transition functions
models{1} = @model_Lopes;
models{2} = @model_Lopes;
models{3} = @model_Lopes;

% Measurement fnctions
Hmodels{1} = [1 0 -1 0 0 0];           % Lopes
Hmodels{2} = [1 0 -1 0 0 0];           % Lopes
Hmodels{3} = [1 0 -1 0 0 0];           % Lopes

% Indices to plot (only potentials, not derivatives)
displayStates{1} = [1 ];
displayStates{2} = [1 ];
displayStates{3} = [1 ];

noiseInd = [2 2 2];     % NEED TO CHECK THESE! (and adjust Wendling noise)
noiseScale = [1 1 1];

crbModels = zeros(6,N,length(models));

priorSigma = [10 100 100];
measureSigma = [10 10 100];

plotStyles={'-','--',':'};

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
P0 = priorSigma(iModel).^2*eye(NStates);

Hfun = @(x)H;
R = measureSigma(iModel)^2*eye(1);
    
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Compute CRLB
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
crbAll = compute_pcrb_P(t,f,F,Hfun,Q,R,m0,P0,100);


crbModels(:,:,iModel) = crbAll;

semilogy(t,sqrt(crbAll(displayInd,:)),plotStyles{iModel})
iModel
hold on;


end

return;
%%
save crb_convergence_results.mat crbModels t priorSigma priorSigma

plotConvergenceResults
