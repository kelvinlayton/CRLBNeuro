
FS = 8;             % fontsize - point
FS_Panel_Label = 12;

font = 'arial';
fig_width = 17;  % cm
fig_height = 7;  % cm

LW = 0.5;         % linewidth

fig_filename = '/Users/klayton/Dropbox/CRLBNeuro/manuscript/figures/pdf/CRBbar.pdf';
figure('units','centimeters','papersize',[fig_height fig_width],'filename',fig_filename,'position',[2 2 fig_width fig_height]);

load crb_ekf_results.mat

%%
mseBar = [];
crbBar = [];
for iModel=1:length(mseModels)
    mseBar = [mseBar; mseModels{iModel}; nan];
    crbBar = [crbBar; crbModels{iModel}; nan];
end


b2=bar(mseBar(:)','BaseValue',1e-6,'BarWidth',0.15); hold on;
b1=bar(crbBar(:)','BaseValue',1e-6); hold on;
set(gca,'YScale','log','YGrid','on')
set(b2,'EdgeColor','none')


ch1 = get(b1,'children');
ch2 = get(b2,'children');
cind = [1;1;1;1;2;2;2;2;3;3;3;3;3;3;4;4;4;4];
set(ch1,'FaceVertexCData',cind)
set(ch2,'FaceVertexCData',cind)
cmap = [0.7556 0.8196 0.4769; 0.6590 0.4332 0.2696; 0.2245 0.3898 0.5888];
cmap = rand(4,3);
colormap(cmap);

for i=1:length(mseBar)
l1=plot(i,mseBar(i),'o');
set(l1,'MarkerSize',7,'MarkerEdgeColor','k','MarkerFaceColor',cmap(cind(i),:));
end

laby=ylabel('MSE');
set(gca,'XTick',[2 6 11 16],'XTickLabel',{'JR','Lopes','Wendling','Linear Sanity Check'})
% set(gca,'YTick',[1e-6 1e-4 1e-2 1e0 1e2])
ylim([1e-6 1e3])

set([gca laby],'fontsize',FS,'fontname',font)

export_fig(fig_filename)