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

% Intial true parameter values
% mode = 'background';
mode = 'alpha';       

NStatesModels = [6 6 10 3];  

% Model parameters
params{1} = SetParametersJR(mode);
params{2} = SetParametersLopes(mode);
params{3} = SetParametersWendling(mode);
params{4} = params{1}; % Only need dt

% Transition functions
models{1} = @model_JR;
models{2} = @model_Lopes;
models{3} = @model_Wendling;
models{4} = @model_Linear;

% Measurement fnctions
Hmodels{1} = [0 0 1 0 -1 0];           % JR
Hmodels{2} = [1 0 -1 0 0 0];           % Lopes
Hmodels{3} = [0 0 1 0 -1 0 -1 0 0 0];  % Wendling
Hmodels{4} = [1 1 1];

% Indices to plot (only potentials, not derivatives)
displayStates{1} = [1 3 5];
displayStates{2} = [1 3 5];
displayStates{3} = [1 3 5 7 9];
displayStates{4} = [1 2 3];
noiseInd = [4 2 4 2];     % NEED TO CHECK THESE! (and adjust Wendling noise)
noiseScale = [1 1 0.4225 1];

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
crbAll = compute_pcrb_P(t,f,F,Hfun,Q,R,m0,P0,100);

crb = crbAll(:,end);


semilogy(t,crbAll(displayInd,:))

% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Compute empirical MSE
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

mse = zeros(NStates,1);
parfor i=1:M
    
    % Generate trajectory realisation
    %
    v = mvnrnd(zeros(NStates,1),Q,N)';

    x = zeros(NStates,N); 
    x(:,1) = mvnrnd(x0,P0)';

    % Euler-Maruyama integration
    %
    for n=1:N-1
        x(:,n+1) = f(x(:,n))  + v(:,n);
    end

    w = mvnrnd(zeros(size(H,1),1),R,N)';
    y = H*x + w;

    % Apply EKF filter
    %
    m = extended_kalman_filter(y,f,F,H,Q,R,m0,P0);

    mse = mse + (m(:,end) - x(:,end)).^2;

end

mse = mse ./ M;

mseModels{iModel} = mse(displayInd);
crbModels{iModel} = crb(displayInd);

end

save crb_ekf_results.mat mseModels crbModels

plotEKFresults
