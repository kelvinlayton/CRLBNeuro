% This code repricates the structure of the Voss / Schiff FN simulation /
% estimation


N = 2000;             	% number of samples
dT = 0.001;          	% sampling time step (global)
dt = 1*dT;            	% integration time step
nn = fix(dT/dt);      	% (used in for loop for forward modelling) the integration time step can be small that the sampling (fix round towards zero)

t = 0:dt:(N-1)*dt;

% Intial true parameter values
mode='alpha';
% mode='seizure';
parameters = SetParametersJR(mode);
parameters.dt = dt;
A=parameters.A;
a=parameters.a;
sigma=parameters.sigma;

% Initialise random number generator for repeatability
%
rng(0);

% Transition model
%
NStates = 6;                           
f = @(x)model_JR(x,'transition',parameters);
F = @(x)model_JR(x,'jacobian',parameters);

% Initialise trajectory state
%
x0 = zeros(NStates,1);                   % initial state
x = zeros(NStates,N); 
x(:,1) = x0;


Q = zeros(NStates);
Q(4,4) = (sqrt(dt)*A*a*sigma)^2;
R = 1^2*eye(1);

% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Generate trajectory
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

% Define input
%
v = mvnrnd(zeros(NStates,1),Q,N)';

% Euler-Maruyama integration
%
for n=1:N-1
    x(:,n+1) = f(x(:,n))  + v(:,n);
end

H = [0 0 1 0 -1 0];           % observation function
w = mvnrnd(zeros(size(H,1),1),R,N)';
y = H*x + w;
figure,plot(t,-H*x)


% save(['JR_' mode '.mat'],'y','t') 


%%
% One column = 17.35 max
% two column = 8.3 cm max
FS = 8;             % fontsize - point
font = 'arial';
fig_width = 5.7;  % cm
fig_height = 2;  % cm
fig_filename = '/Users/dean/Projects/CRLBNeuro/manuscript/figures/pdf/JR_seizure.pdf';
figure('units','centimeters','papersize',[fig_height fig_width],'filename',fig_filename,'position',[2 2 fig_width fig_height]);

LW = 1;         % linewidth
plot(t,-y,'k','linewidth',LW)
axis tight
axis off
set(gca,'fontsize',FS,'fontname',font)
text_x = -0.25;
text_y = max(-y);
text(text_x,text_y,'D','fontsize',FS,'fontname',font)
% xlabel('Time (s)','fontsize',FS,'fontname',font)
% ylabel('Amplitude','fontsize',FS,'fontname',font)

%% 
% % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% % Run EKF
% % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% 
% % Prior distribution (defined by m0 & P0)
% %
% m0 = x0;
% P0 = 100.^2*eye(NStates);
% 
% 
% % Apply EKF filter
% %
% m = extended_kalman_filter(y,f,F,H,Q,R,m0,P0);
% 
% figure
% ax1=subplot(211);
% plot(t,x([3 5],:)'); hold on;
% plot(t,m([3 5],:)','--');
% ax2=subplot(212);
% plot(t,y)
% 
% linkaxes([ax1 ax2],'x');
% 
% return
% 
% %% 
% % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% % Compute the posterior Cramer-Rao bound (PCRB)
% % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% 
% M = 100;    % Number of Monte Carlo samples
% 
% pcrb = compute_pcrb_P(t,f,F,@(x)H,Q,R,m0,P0,M);
% 
% 
% %% 
% % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% % Compute the MSE of the extended Kalman filter 
% % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% 
% %
% num_trials = 10;
% error = zeros(NStates,N);
% 
% parfor r=1:num_trials
%     
%     % Create new trajectory realisation
%     %
%     v = mvnrnd(zeros(NStates,1),Q,N)';
%     x = zeros(NStates,N);
%     x(:,1)=mvnrnd(m0,P0)';
%     for i=2:N
%         x(:,i) = f(x(:,i-1)) +  v(:,i);
%     end
% 
%     % Generate new observations 
%     %
%     w = mvnrnd(zeros(1,size(H,1)),R,N)';
%     z = H*x + w;
% 
%     % Apply EKF filter
%     %
%     m = extended_kalman_filter(z,f,F,H,Q,R,m0,P0);
% 
%     % Accumulate the estimation error
%     %
%     error = error + (x-m).^2;
% end
% 
% % Calculate the mean squared error
% %
% mse = error ./ num_trials;
% 
% %% Plot MSE and the PCRB
% %
% figure
% % semilogy(t,sqrt(sum(mse([3 5],:))),'x-')
% % hold on;
% semilogy(t,sqrt(pcrb([3 5],:)),'-','LineWidth',2);
% title(mode);
% 
% grid on;
% xlabel('Time (s)');
% ylabel('RMSE');
% % export_fig(['../figures/nmm_pcrb_' mode '.pdf'])
