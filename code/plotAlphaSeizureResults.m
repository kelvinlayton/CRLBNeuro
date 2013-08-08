

FS = 8;             % fontsize - point
FS_Panel_Label = 12;

font = 'arial';
fig_width = 17;  % cm
fig_height = 7;  % cm

LW = 0.5;         % linewidth

fig_filename = '/Users/kelvin/Dropbox/CRLBNeuro/manuscript/figures/pdf/AlphaSeizureBar.pdf';
figure('units','centimeters','papersize',[fig_height fig_width],'filename',fig_filename,'position',[2 2 fig_width fig_height]);

load crb_seizure_results.mat

%%
crbBar = [];
for iModel=[1 3 2 4]
    crbBar = [crbBar; sum(sqrt(crbModels{iModel})); nan];
end

b1=bar(crbBar(:)','BaseValue',1e-6); hold on;
% set(gca,'YScale','log','YGrid','on')
set(gca,'YScale','linear','YGrid','on')
% xlim([0 20])

ch1 = get(b1,'children');
% cind = [1;1;1;1;2;2;2;2;2;2;3;3;3;3;4;4;4;4;4;4];
cind = [1 1 2 2 3 3 4 4]';
set(ch1,'FaceVertexCData',cind)
cmap = [0.7556 0.8196 0.4769; 0.6590 0.4332 0.2696; 0.2245 0.3898 0.5888; 0.1831 0.9336 0.5528];
% cmap = rand(4,3);
colormap(cmap);

laby=ylabel('RMSE (mV)');
% set(gca,'XTick',[2 6 11 16],'XTickLabel',{'JR alpha','JR seizure','Wendling alpha','Wendling seizure'})
set(gca,'XTick',[1 3 5 7],'XTickLabel',{'JR alpha','JR seizure','Wendling alpha','Wendling seizure'})
% ylim([1e-6 1e3])
ylim([0 10])

set([gca laby],'fontsize',FS,'fontname',font)

return
%%
export_fig(fig_filename)