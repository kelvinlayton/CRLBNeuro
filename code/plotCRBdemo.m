
FS = 8;             % fontsize - point
FS_Panel_Label = 12;

font = 'arial';
fig_width = 20;  % cm
fig_height = 20;  % cm

LW = 0.5;         % linewidth

fig_filename = '/Users/klayton/Dropbox/CRLBNeuro/manuscript/figures/pdf/CRBbar.pdf';
figure('units','centimeters','papersize',[fig_height fig_width],'filename',fig_filename,'position',[2 2 fig_width fig_height]);

xValues=linspace(-2,8,100);
yValues=linspace(-4,6,100);
[X,Y]=meshgrid(xValues,yValues);

x0 = 3;
R = 1;
P0 = 4;

MU = X./100;
H = 1./100.*ones(size(MU));


p_like = normpdf(Y,MU,R);
p_x = normpdf(X,x0,P0);
p_xy = p_like.*p_x;

der_p_x = -inv(P0)*ones(size(X));
der_p_like = -H.^2*inv(R);
der_p_xy = der_p_x + der_p_like;
% der_p_xy = der_p_xy.*p_xy;

R = 0.5;
P0 = 1;

p_like = normpdf(Y,MU,R);
p_x = normpdf(X,x0,P0);
p_xy2 = p_like.*p_x;

der_p_x = -inv(P0)*ones(size(X));
der_p_like = -H.^2*inv(R);
der_p_xy2 = der_p_x + der_p_like;
% der_p_xy2 = der_p_xy2.*p_xy2;

axes('Pos',[0 0 0.5 0.5])
mesh(xValues,yValues,p_xy)

axes('Pos',[0.5 0 0.5 0.5])
mesh(xValues,yValues,-der_p_xy)


axes('Pos',[0 0.5 0.5 0.5])
mesh(xValues,yValues,p_xy2)

axes('Pos',[0.5 0.5 0.5 0.5])
mesh(xValues,yValues,-der_p_xy2)

return

%%
mseBar = [];
crbBar = [];
for iModel=1:length(mseModels)
    mseBar = [mseBar; sqrt(mseModels{iModel}); nan];
    crbBar = [crbBar; sqrt(crbModels{iModel}); nan];
end


b2=bar(mseBar(:)','BaseValue',1e-6,'BarWidth',0.15); hold on;
b1=bar(crbBar(:)','BaseValue',1e-6); hold on;
set(gca,'YScale','log','YGrid','on')
% set(gca,'YScale','linear','YGrid','on')
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

laby=ylabel('RMSE (mV)');
set(gca,'XTick',[2 6 11 16],'XTickLabel',{'JR','Lopes','Wendling','Linear Sanity Check'})
% set(gca,'YTick',[1e-6 1e-4 1e-2 1e0 1e2])
ylim([1e-6 1e3])
% ylim([0 10])

set([gca laby],'fontsize',FS,'fontname',font)

return
%%
export_fig(fig_filename)