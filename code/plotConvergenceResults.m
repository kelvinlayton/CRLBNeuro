

FS = 8;             % fontsize - point
FS_Panel_Label = 12;

font = 'arial';
fig_width = 17;  % cm
fig_height = 7;  % cm

LW = 0.5;         % linewidth

fig_filename = '/Users/klayton/Dropbox/CRLBNeuro/manuscript/figures/pdf/DynamicCRB.pdf';
figure('units','centimeters','papersize',[fig_height fig_width],'filename',fig_filename,'position',[2 2 fig_width fig_height]);

load crb_convergence_results.mat

%%
displayInd = 1;

crbData = squeeze(sqrt(crbModels(displayInd,2:end,:)));
lHandle=plot(1e3*t(2:end),crbData);

set(lHandle,'LineWidth',2)
set(lHandle(1),'Color',[60 117 184]./255,'LineStyle','-');
set(lHandle(2),'Color',[239 133 0]./255,'LineStyle','--');
set(lHandle(3),'Color',[198 54 34]./255,'LineStyle','-.');
set(gca,'Position',[0.1 0.17 0.8 0.8],'TickDir','out');
xlim([0 80]);
% grid on;
box off

labx=xlabel('Time (ms)');
laby=ylabel('RMSE (mV)');
legText={};
for i=1:length(priorSigma)
legText = cat(1,legText,sprintf('$P_0=%d, R=%d$',priorSigma(i)^2,measureSigma(i)^2));
end
leg=legend(legText);

set([gca laby labx leg],'fontsize',FS,'fontname',font)
set(leg,'interpreter','latex','Position',get(leg,'Position').*[0.95 1 1.1 1.1]);
legend('boxoff')
return
%%
export_fig(fig_filename)