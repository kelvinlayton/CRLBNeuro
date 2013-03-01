
% This script performs a quick visual test of the calculated jacobian
% 
% A single parameter is changed and the output vector is plotted as well as
% the linear approximation given by the jacobian
%

dt = 0.001; 
params = SetParametersJR('alpha');
params.dt = dt;

NStates = 6;                           
f = @(x)model_JR_erf(x,'transition',params);
F = @(x)model_JR_erf(x,'jacobian',params);

rng(0);

p=linspace(-1,1,1000);
P=length(p);

x=repmat(rand(6,1),1,P);
y=zeros(NStates,P);
yhat=zeros(NStates,P);

% Calculate the output vector as the parameter is varied
%
perturbIdx=1;
x(perturbIdx,:)=p;
for i=1:P
    y(:,i)=f(x(:,i));
end

% Calculate the linear approximation using the jacobian about the specified
% fixed point
%
p0=0.06;
[~,idx] = min(abs(p-p0));
p0=p(idx);
x0=x(:,idx);
for i=1:P
    yhat(:,i)=f(x0) + F(x0)*(x(:,i)-x0) ;
end

% Plot elements of vector for different paramter values
%
plot(p,y);
hold on;
plot(p0,y(:,idx),'o')
plot(p,yhat,'--')
xlim(0.2*[-1 1]);
