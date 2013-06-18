clc
close all
clear

t=linspace(0,0.1,1000);
A = 3.25;       % mV
a = 100;        % s^{-1}
h_e=A*a*t.*exp(-a*t);

fig_height = 3;
fig_width = 3;
FS_label = 8;
FS_ticks = 8;

set(0,'defaultAxesFontName', 'timesnewroman')
set(0,'defaultTextFontName', 'timesnewroman')
set(0,'defaultUicontrolFontName', 'timesnewroman')
set(0,'defaultUitableFontName', 'timesnewroman')
set(0,'defaultUipanelFontName', 'timesnewroman')
set(0,'defaultuicontrolunits','centimeters') ;

%%

figure('units','centimeters',...
    'color','white',...
    'position',[0 0 fig_width fig_height],...
    'papersize',[fig_width fig_height],...
    'PaperPositionMode','auto',...
    'renderer','painters')

ax1 = axes;
Left =1.5;
Bottom = 0.9;
Width = 2;
Height = 2;
set(ax1,'units','centimeters','position',[Left Bottom Width Height]);

plot(ax1,[1/a 1/a],[0 max(h_e)],'-.r'),hold on
plot(ax1,[t(1) t(end)],[max(h_e) max(h_e)],'--')
plot(ax1,t,h_e,'k')
% axis tight
% plot(t,(A/a)*ones(1,1000),'k')

box off
set(gca,'fontsize',FS_ticks,'xticklabel',{},'yticklabel',{})
leg = legend('$\tau$','$\propto \alpha$');
set(leg,'units','centimeters',...
    'interpreter','latex',...
    'location','southeast',...
    'box','off',...
    'fontsize',FS_label,...
    'fontname','timesnewroman',...
    'position',[2.2 1.5 1.25 1]);
axis tight
ylim([0 max(h_e)+0.1])
xlabel('Time (s)','fontname','arial','fontsize',FS_label)
ylabel('PSRK (mV)','fontsize',FS_label,'fontname','arial')



%%

e_0 = 1;          % s^{-1}
v_0 = 6;              % mV
r = 0.56;           % mV^{-1}

min_v = 15;
max_v = 20;
v = linspace(v_0-min_v,v_0+20,1000);

f_v = 1 ./ (1 + exp(r*(v_0 - v)));    

% f_v = f(v);
figure('units','centimeters',...
    'color','white',...
    'position',[2 2 fig_width fig_height],...
    'papersize',[fig_width fig_height],...
    'PaperPositionMode','auto',...
    'renderer','painters')
ax2 = axes;
Left =1.5;
Bottom = 0.9;
Width = 2;
Height = 2;
set(ax2,'units','centimeters','position',[Left Bottom Width Height]);

plot(ax2,[v_0 v_0],[-0.05 1.05],'r--'),hold on
plot(ax2,v,f_v,'k')
xlabel('Membrane Pot. (mV)','fontsize',FS_label,'fontname','arial')
ylabel('Rate (APs/s)','fontsize',FS_label,'fontname','arial')
box off
ylim([-0.25 5.25])
xlim([v_0-min_v,v_0+20])
set(gca,'fontsize',FS_ticks,'fontname','arial','xticklabel',{},'yticklabel',{})
% plot(v,2*e_0*ones(1,1000),'r')
leg = legend('$v_0$');
set(leg,'units','centimeters',...
    'interpreter','latex',...
    'location','southeast',...
    'box','off',...
    'fontsize',FS_label,...
    'fontname','arial',...
    'position',[2.8 1.25 0.8 0.1]);
axis tight
