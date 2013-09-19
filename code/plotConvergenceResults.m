

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
l1=plot(t(2:end),crbData,'-');

% set(gca,'YScale','log','YGrid','on')
% set(gca,'YScale','linear','YGrid','on')
% xlim([0 20])
grid on;


laby=ylabel('RMSE (mV)');
% set(gca,'XTick',[2 6 11 16],'XTickLabel',{'JR alpha','JR seizure','Wendling alpha','Wendling seizure'})
% set(gca,'XTick',[1 3 5 7],'XTickLabel',{'JR alpha','JR seizure','Wendling alpha','Wendling seizure'})
% ylim([1e-6 1e3])
% ylim([0 10])

set([gca laby],'fontsize',FS,'fontname',font)

return
%%
export_fig(fig_filename)