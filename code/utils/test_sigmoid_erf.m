
% This script determines the relationship between the width of the 
% logistic sigmoid function and the 'equivalent' width of a sigmoid 
% defined using the error function.
% 

x = linspace(-10,10,1000);

v0Sig = 0;
v0Erf = 0;
rSig=1;
rErf = find_erf_width(rSig);

sig1 = 1 ./ (1 + exp(rSig*(v0Sig - x)));
sig2 = 0.5 .* ( 1 + erf((x - v0Erf)./(sqrt(2)*rErf)) );

figure
plot(x,sig1,'r');
hold on;
plot(x,sig2,'b');
legend('Sigmoid','Error function');

%% Determine relationship between width of sigmoid functions

rSigValues = linspace(0.01,20,100)';

rErfValues = zeros(size(rSigValues));
for iR=1:length(rSigValues)
    rSig = rSigValues(iR);
    
    rErfValues(iR) = find_erf_width(rSig);

end

% Linear fit to inverse
%
fitType=fittype('a*x');
linearFit = fit(rSigValues,1./rErfValues,fitType)

% Display results
%
fprintf('=======================================\n');
fprintf('(width erf) = %.3f / (width logistic)\n',1./linearFit.a)
fprintf('=======================================\n');

% Plot
%
figure
plot(rSigValues,1./rErfValues,'-x')
hold on;
plot(rSigValues,(linearFit.a .* rSigValues),'r');
xlabel('Logistic width'); ylabel('1 / (ERF width)');

% plot(rSigValues,rErfValues,'-x')
% hold on;
% plot(rSigValues,1./(linearFit.a .* rSigValues),'r');
% xlabel('Sigmoid width'); ylabel('ERF width');
