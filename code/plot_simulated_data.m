% One column = 17.35 max
% two column = 8.3 cm max
close all
clc
clear

FS = 8;             % fontsize - point
FS_Panel_Label = 12;

font = 'arial';
fig_width = 17;  % cm
fig_height = 5;  % cm

LW = 0.5;         % linewidth

min_y = -18;
max_y = 27;

fig_filename = '/Users/dean/Projects/CRLBNeuro/manuscript/figures/pdf/JR_seizure.pdf';
figure('units','centimeters','papersize',[fig_height fig_width],'filename',fig_filename,'position',[2 2 fig_width fig_height]);


%
load Lopes_alpha.mat
y = detrend(y);
max_l_a = max(y(3001:end))
min_l_a = min(y(3001:end))

ax(1) = subplot(231);

plot(t(3001:end)-3,y(3001:end),'k','linewidth',LW)

ylim(ax(1),[min_y,max_y])
xlim(ax(1),[0 2])

axis off
set(gca,'fontsize',FS,'fontname',font)
text_x = -0.25;
text_y = max(y(3001:end));
text(text_x,text_y,'A','fontsize',FS,'fontname',font)

%
load JR_alpha.mat
y = detrend(y);
max_j_a = max(y)
min_j_a = min(y)

ax(2) = subplot(232);

plot(t,y,'k','linewidth',LW)

ylim(ax(2),[min_y,max_y])
xlim(ax(2),[0 2])

axis off
set(gca,'fontsize',FS,'fontname',font)
text_x = -0.25;
text_y = max(y);
text(text_x,text_y,'B','fontsize',FS,'fontname',font)

load JR_seizure.mat
y = detrend(y);
max_j_s = max(y)
min_j_s = min(y)

ax(3) = subplot(235);

plot(t,y,'k','linewidth',LW)

ylim(ax(3),[min_y,max_y])
xlim(ax(3),[0 2])

axis off
set(gca,'fontsize',FS,'fontname',font)
text_x = -0.25;
text_y = max(-y);
text(text_x,text_y,'D','fontsize',FS,'fontname',font)

%
load Wendling_alpha.mat
y = detrend(y);
max_w_a = max(y)
min_w_a = min(y)

ax(4) = subplot(233);

plot(t,y,'k','linewidth',LW)

ylim(ax(4),[min_y,max_y])
xlim(ax(4),[0 2])

axis off
set(gca,'fontsize',FS,'fontname',font)
text_x = -0.25;
text_y = max(y);
text(text_x,text_y,'C','fontsize',FS,'fontname',font)

load Wendling_spikes.mat
y = detrend(y);
max_w_s = max(y)
min_w_s = min(y)

ax(5) = subplot(236);

plot(t,y,'k','linewidth',LW)

ylim(ax(5),[min_y,max_y])
xlim(ax(5),[0 2])

axis off
set(gca,'fontsize',FS,'fontname',font)
text_x = -0.25;
text_y = max(y);
text(text_x,text_y,'E','fontsize',FS_Panel_Label,'fontname',font)

% this is the legenend for time and amplitude
ax(6) = subplot(234);

ylim(ax(6),[min_y,max_y])
xlim(ax(6),[0 2])

axis off
line([0 1],[0 0],'color','black')
line([0 0],[0 20],'color','black')
text(.4,-max_y/4,'1 s','fontsize',FS,'fontname',font)
text(-.15,0,'20 mV','rotation',90,'fontsize',FS,'fontname',font)
% text(-.15,-max_y/4,'Amplitude','rotation',90)