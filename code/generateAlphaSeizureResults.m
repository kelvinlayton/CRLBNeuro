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


NStatesModels = [6 10 6 10];  

% Model parameters
params{1} = SetParametersJR('alpha');
params{2} = SetParametersWendling('alpha');
params{3} = SetParametersJR('seizure');
params{4} = SetParametersWendling('spikes');

% Transition functions
models{1} = @model_JR;
models{2} = @model_Wendling;
models{3} = @model_JR;
models{4} = @model_Wendling;

% Measurement fnctions
Hmodels{1} = [0 0 1 0 -1 0];           % JR
Hmodels{2} = [0 0 1 0 -1 0 -1 0 0 0];  % Wendling
Hmodels{3} = [0 0 1 0 -1 0];           % JR
Hmodels{4} = [0 0 1 0 -1 0 -1 0 0 0];  % Wendling

% Indices to plot (only potentials, not derivatives)
displayStates{1} = [1 3 5];
displayStates{2} = [1 3 5 7 9];
displayStates{3} = [1 3 5];
displayStates{4} = [1 3 5 7 9];

noiseInd = [4 4 4 4 ];     % NEED TO CHECK THESE! (and adjust Wendling noise)
noiseScale = [1 0.4225 1 0.4225];

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
P0 = 10.^2*eye(NStates);


Hfun = @(x)H;
R = 10^2*eye(1);
    
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Compute CRLB
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
crbAll = compute_pcrb_P(t,f,F,Hfun,Q,R,m0,P0,500);

crb = crbAll(:,end);


semilogy(t,crbAll(displayInd,:))

crbModels{iModel} = crb(displayInd);

end

save crb_seizure_results.mat crbModels

plotAlphaSeizureResults
