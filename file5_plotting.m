
%% Clear up

close all
clear
clc

%% Add Path

addpath(genpath('/Users/danielalexander/Documents/Analysis/DQE/DQE Analysis'));
addpath(genpath('/Users/danielalexander/Documents/Matlab Functions'));

%% Initialize

load('/Users/danielalexander/Documents/Analysis/DQE/DQE Analysis/Results/DQE_all_results.mat');

colors = linspecer(8);

% cmos_gains = {'100','176','314','418','554','740','971'};
cmos_gains = {'8%', '14%','25%','32%','43%','57%','76%','100%'};

% red_int_gains = {'100','129','162','226','317'};
red_int_gains = {'32%','41%','51%','71%','100%'};
% LN_int_gains = {'100','145','197','275','347'};
LN_int_gains = {'29%','42%','57%','79%','100%'};

l = 80;
pixelsize = 30/l;

Fs = 1/pixelsize;
N = length(r)*2; % to account for positive & negative
freqdom = (Fs/2*0:Fs/N:Fs/2-Fs/N)';

ms = 6; % markersize
lw = 1; % Linewidth
start_index = 2; % might want to avoid first point?
xmin = 0; xmax = 1;


figSize1 = [100,100,1000,420];
figSize2 = [500,500,500,420];

labelPos = [-.2, 1.05];

%% MTF

figure('position', figSize2);

hold on;
plot(MTF.freq, MTF.RedCam.nofilt.mtf,   '^-',  'Color', colors(1,:), 'LineWidth', lw, 'MarkerSize', ms);
plot(MTF.freq, MTF.RedCam.filt.mtf,     '^--', 'Color', colors(1,:), 'LineWidth', lw, 'MarkerSize', ms, 'MarkerFaceColor', colors(1,:));
plot(MTF.freq, MTF.LNCam.nofilt.mtf,   's-',  'Color', colors(2,:), 'LineWidth', lw, 'MarkerSize', ms);
plot(MTF.freq, MTF.LNCam.filt.mtf,     's--', 'Color', colors(2,:), 'LineWidth', lw, 'MarkerSize', ms, 'MarkerFaceColor', colors(2,:));
yline(0.1, 'k:', 'Linewidth', 1);
ylabel('MTF');
xlabel('Frequency (mm^{-1})');
xlim([xmin, xmax]);
legend('Gen3 Raw','Gen3 Filtered', 'Gen2+ Raw','Gen2+ Filtered');
legend boxoff
grid on;
ax = gca;
ax.FontSize = 16; 

cd '/Users/danielalexander/Documents/Analysis/DQE/DQE Analysis/'
saveas(gcf, 'Figs/MTF.png')

% Calculate f_10

f10_redNoFilt   = interp1(MTF.RedCam.nofilt.mtf,  MTF.freq, 0.1);
disp(['f10_redNoFilt = ', num2str(f10_redNoFilt)]);
f10_redFilt     = interp1(MTF.RedCam.filt.mtf,    MTF.freq, 0.1);
disp(['f10_redFilt = ', num2str(f10_redFilt)]);

f10_lnNoFilt    = interp1(MTF.LNCam.nofilt.mtf,   MTF.freq,  0.1);
disp(['f10_lnNoFilt = ', num2str(f10_lnNoFilt)]);
f10_lnFilt      = interp1(MTF.LNCam.filt.mtf,     MTF.freq,  0.1);
disp(['f10_lnFilt = ', num2str(f10_lnFilt)]);

%% NPS - Vary CMOS Gain

figure('position', figSize1);

subplot(1,2,1);
text(labelPos(1),labelPos(2),'a)', 'Units', 'normalized', 'Fontsize', 16)
hold on;
plot(freqdom(start_index:end), NPS1D_RedCam_avg.i1p602.a425(start_index:end), '^-', 'Color', colors(1,:), 'LineWidth', lw, 'MarkerSize', ms);
plot(freqdom(start_index:end), NPS1D_RedCam_avg.i1p602.a400(start_index:end), 'o-', 'Color', colors(2,:), 'LineWidth', lw, 'MarkerSize', ms);
plot(freqdom(start_index:end), NPS1D_RedCam_avg.i1p602.a375(start_index:end), 'd-', 'Color', colors(3,:), 'LineWidth', lw, 'MarkerSize', ms);
plot(freqdom(start_index:end), NPS1D_RedCam_avg.i1p602.a350(start_index:end), 'v-', 'Color', colors(4,:), 'LineWidth', lw, 'MarkerSize', ms);
plot(freqdom(start_index:end), NPS1D_RedCam_avg.i1p602.a325(start_index:end), 's-', 'Color', colors(5,:), 'LineWidth', lw, 'MarkerSize', ms);
plot(freqdom(start_index:end), NPS1D_RedCam_avg.i1p602.a300(start_index:end), '+-', 'Color', colors(6,:), 'LineWidth', lw, 'MarkerSize', ms);
plot(freqdom(start_index:end), NPS1D_RedCam_avg.i1p602.a250(start_index:end), '*-', 'Color', colors(7,:), 'LineWidth', lw, 'MarkerSize', ms);
plot(freqdom(start_index:end), NPS1D_RedCam_avg.i1p602.a200(start_index:end), 'p-', 'Color', colors(8,:), 'LineWidth', lw, 'MarkerSize', ms);
legend(flip(cmos_gains), 'FontSize', 16, 'NumColumns',2);
legend boxoff
ylabel('NPS (mm^2)');
xlabel('Frequency (mm^{-1})');
xlim([xmin,xmax]);
ylim([1e2,2e6]);
ax = gca;
ax.FontSize = 16; 
set(gca, 'YScale', 'log')
grid on;

subplot(1,2,2);
text(labelPos(1),labelPos(2),'b)', 'Units', 'normalized', 'Fontsize', 16)
hold on;
plot(freqdom(start_index:end), NPS1D_LNCam_avg.i1p848.a425(start_index:end), '^-', 'Color', colors(1,:), 'LineWidth', lw, 'MarkerSize', ms);
plot(freqdom(start_index:end), NPS1D_LNCam_avg.i1p848.a400(start_index:end), 'o-', 'Color', colors(2,:), 'LineWidth', lw, 'MarkerSize', ms);
plot(freqdom(start_index:end), NPS1D_LNCam_avg.i1p848.a375(start_index:end), 'd-', 'Color', colors(3,:), 'LineWidth', lw, 'MarkerSize', ms);
plot(freqdom(start_index:end), NPS1D_LNCam_avg.i1p848.a350(start_index:end), 'v-', 'Color', colors(4,:), 'LineWidth', lw, 'MarkerSize', ms);
plot(freqdom(start_index:end), NPS1D_LNCam_avg.i1p848.a325(start_index:end), 's-', 'Color', colors(5,:), 'LineWidth', lw, 'MarkerSize', ms);
plot(freqdom(start_index:end), NPS1D_LNCam_avg.i1p848.a300(start_index:end), '+-', 'Color', colors(6,:), 'LineWidth', lw, 'MarkerSize', ms);
plot(freqdom(start_index:end), NPS1D_LNCam_avg.i1p848.a250(start_index:end), '*-', 'Color', colors(7,:), 'LineWidth', lw, 'MarkerSize', ms);
plot(freqdom(start_index:end), NPS1D_LNCam_avg.i1p848.a200(start_index:end), 'p-', 'Color', colors(8,:), 'LineWidth', lw, 'MarkerSize', ms);
legend(flip(cmos_gains), 'FontSize', 16, 'NumColumns',2);
legend boxoff
ylabel('NPS (mm^2)');
xlabel('Frequency (mm^{-1})');
xlim([xmin,xmax]);
ylim([1e2,2e6]);
ax = gca;
ax.FontSize = 16; 
set(gca, 'YScale', 'log')
grid on;

cd '/Users/danielalexander/Documents/Analysis/DQE/DQE Analysis/'
saveas(gcf, 'Figs/NPS_CMOS.png')

%% NPS - Vary Intensifier Gain

figure('position', figSize1);

subplot(1,2,1);
text(labelPos(1),labelPos(2),'a)', 'Units', 'normalized', 'Fontsize', 16)
hold on;
plot(freqdom(start_index:end), NPS1D_RedCam_avg.i1p602.a425(start_index:end), '^-', 'Color', colors(1,:), 'LineWidth', lw, 'MarkerSize', ms);
plot(freqdom(start_index:end), NPS1D_RedCam_avg.i1p717.a425(start_index:end), 'o-', 'Color', colors(2,:), 'LineWidth', lw, 'MarkerSize', ms);
plot(freqdom(start_index:end), NPS1D_RedCam_avg.i1p847.a425(start_index:end), 'd-', 'Color', colors(3,:), 'LineWidth', lw, 'MarkerSize', ms);
plot(freqdom(start_index:end), NPS1D_RedCam_avg.i1p923.a425(start_index:end), 'v-', 'Color', colors(4,:), 'LineWidth', lw, 'MarkerSize', ms);
plot(freqdom(start_index:end), NPS1D_RedCam_avg.i2p016.a425(start_index:end), 's-', 'Color', colors(5,:), 'LineWidth', lw, 'MarkerSize', ms);
legend(flip(red_int_gains), 'FontSize', 16);
legend boxoff
ylabel('NPS (mm^2)');
xlabel('Frequency (mm^{-1})');
xlim([xmin,xmax]);
ylim([4e3,2e6]);
ax = gca;
ax.FontSize = 16; 
set(gca, 'YScale', 'log')
grid on;

subplot(1,2,2);
text(labelPos(1),labelPos(2),'b)', 'Units', 'normalized', 'Fontsize', 16)
hold on;
plot(freqdom(start_index:end), NPS1D_LNCam_avg.i1p848.a425(start_index:end), '^-', 'Color', colors(1,:), 'LineWidth', lw, 'MarkerSize', ms);
plot(freqdom(start_index:end), NPS1D_LNCam_avg.i1p922.a425(start_index:end), 'o-', 'Color', colors(2,:), 'LineWidth', lw, 'MarkerSize', ms);
plot(freqdom(start_index:end), NPS1D_LNCam_avg.i2p034.a425(start_index:end), 'd-', 'Color', colors(3,:), 'LineWidth', lw, 'MarkerSize', ms);
plot(freqdom(start_index:end), NPS1D_LNCam_avg.i2p135.a425(start_index:end), 'v-', 'Color', colors(4,:), 'LineWidth', lw, 'MarkerSize', ms);
plot(freqdom(start_index:end), NPS1D_LNCam_avg.i2p256.a425(start_index:end), 's-', 'Color', colors(5,:), 'LineWidth', lw, 'MarkerSize', ms);
legend(flip(LN_int_gains), 'FontSize', 16);
legend boxoff
ylabel('NPS (mm^2)');
xlabel('Frequency (mm^{-1})');
xlim([xmin,xmax]);
ylim([4e3,2e6]);
ax = gca;
ax.FontSize = 16; 
set(gca, 'YScale', 'log')
grid on;

cd '/Users/danielalexander/Documents/Analysis/DQE/DQE Analysis/'
saveas(gcf, 'Figs/NPS_Int.png')

%% NPS - Compare Filtered/Unfiltered

figure('position', figSize2);

hold on;
plot(freqdom(start_index:end), NPS1D_RedCam_avg.i1p602.a425(start_index:end),       '^-',  'Color', colors(1,:), 'LineWidth', lw, 'MarkerSize', ms, 'MarkerFaceColor', colors(1,:));
plot(freqdom(start_index:end), NPS1D_RedCam_filt_avg.i1p602.a425(start_index:end),  '^--', 'Color', colors(1,:), 'LineWidth', lw, 'MarkerSize', ms);
plot(freqdom(start_index:end), NPS1D_LNCam_avg.i1p848.a425(start_index:end),        's-',  'Color', colors(2,:), 'LineWidth', lw, 'MarkerSize', ms, 'MarkerFaceColor', colors(2,:));
plot(freqdom(start_index:end), NPS1D_LNCam_filt_avg.i1p848.a425(start_index:end),   's--', 'Color', colors(2,:), 'LineWidth', lw, 'MarkerSize', ms);
set(gca, 'YScale', 'log')
ylabel('NPS (mm^2)');
xlabel('Frequency (mm^{-1})');
xlim([xmin,xmax]);
legend('Gen3 Raw', 'Gen3 Filtered', 'Gen2+ Raw', 'Gen2+ Filtered');
legend boxoff
grid on;
ax = gca;
ax.FontSize = 16; 

cd '/Users/danielalexander/Documents/Analysis/DQE/DQE Analysis/'
saveas(gcf, 'Figs/NPS_Comp.png')


%% DQE - Vary CMOS Gain

figure('position', figSize1);

subplot(1,2,1);
text(labelPos(1),labelPos(2),'a)', 'Units', 'normalized', 'Fontsize', 16)
hold on;
plot(freqdom(start_index:end), 100*DQE_RedCam.i1p602.a425(start_index:end), '^-', 'Color', colors(1,:), 'LineWidth', lw, 'MarkerSize', ms);
plot(freqdom(start_index:end), 100*DQE_RedCam.i1p602.a400(start_index:end), 'o-', 'Color', colors(2,:), 'LineWidth', lw, 'MarkerSize', ms);
plot(freqdom(start_index:end), 100*DQE_RedCam.i1p602.a375(start_index:end), 'd-', 'Color', colors(3,:), 'LineWidth', lw, 'MarkerSize', ms);
plot(freqdom(start_index:end), 100*DQE_RedCam.i1p602.a350(start_index:end), 'v-', 'Color', colors(4,:), 'LineWidth', lw, 'MarkerSize', ms);
plot(freqdom(start_index:end), 100*DQE_RedCam.i1p602.a325(start_index:end), 's-', 'Color', colors(5,:), 'LineWidth', lw, 'MarkerSize', ms);
plot(freqdom(start_index:end), 100*DQE_RedCam.i1p602.a300(start_index:end), '+-', 'Color', colors(6,:), 'LineWidth', lw, 'MarkerSize', ms);
plot(freqdom(start_index:end), 100*DQE_RedCam.i1p602.a250(start_index:end), '*-', 'Color', colors(7,:), 'LineWidth', lw, 'MarkerSize', ms);
plot(freqdom(start_index:end), 100*DQE_RedCam.i1p602.a200(start_index:end), 'p-', 'Color', colors(8,:), 'LineWidth', lw, 'MarkerSize', ms);
legend(flip(cmos_gains), 'FontSize', 16, 'NumColumns',1);
legend boxoff
ylabel('DQE (%)');
xlabel('Frequency (mm^{-1})');
xlim([xmin,xmax]);
ylim([0,1.7]);
ax = gca;
ax.FontSize = 16; 
grid on;

subplot(1,2,2);
text(labelPos(1),labelPos(2),'b)', 'Units', 'normalized', 'Fontsize', 16)
hold on;
plot(freqdom(start_index:end), 100*DQE_LNCam.i1p848.a425(start_index:end), '^-', 'Color', colors(1,:), 'LineWidth', lw, 'MarkerSize', ms);
plot(freqdom(start_index:end), 100*DQE_LNCam.i1p848.a400(start_index:end), 'o-', 'Color', colors(2,:), 'LineWidth', lw, 'MarkerSize', ms);
plot(freqdom(start_index:end), 100*DQE_LNCam.i1p848.a375(start_index:end), 'd-', 'Color', colors(3,:), 'LineWidth', lw, 'MarkerSize', ms);
plot(freqdom(start_index:end), 100*DQE_LNCam.i1p848.a350(start_index:end), 'v-', 'Color', colors(4,:), 'LineWidth', lw, 'MarkerSize', ms);
plot(freqdom(start_index:end), 100*DQE_LNCam.i1p848.a325(start_index:end), 's-', 'Color', colors(5,:), 'LineWidth', lw, 'MarkerSize', ms);
plot(freqdom(start_index:end), 100*DQE_LNCam.i1p848.a300(start_index:end), '+-', 'Color', colors(6,:), 'LineWidth', lw, 'MarkerSize', ms);
plot(freqdom(start_index:end), 100*DQE_LNCam.i1p848.a250(start_index:end), '*-', 'Color', colors(7,:), 'LineWidth', lw, 'MarkerSize', ms);
plot(freqdom(start_index:end), 100*DQE_LNCam.i1p848.a200(start_index:end), 'p-', 'Color', colors(8,:), 'LineWidth', lw, 'MarkerSize', ms);
legend(flip(cmos_gains), 'FontSize', 16, 'NumColumns',1);
legend boxoff
ylabel('DQE (%)');
xlabel('Frequency (mm^{-1})');
xlim([xmin,xmax]);
ylim([0,1.7]);
ax = gca;
ax.FontSize = 16; 
grid on;

cd '/Users/danielalexander/Documents/Analysis/DQE/DQE Analysis/'
saveas(gcf, 'Figs/DQE_CMOS.png')


%% DQE - Vary Intensifier Gain

figure('position', figSize1);

subplot(1,2,1);
text(labelPos(1),labelPos(2),'a)', 'Units', 'normalized', 'Fontsize', 16)
hold on;
plot(freqdom(start_index:end), 100*DQE_RedCam.i1p602.a425(start_index:end), '^-', 'Color', colors(1,:), 'LineWidth', lw, 'MarkerSize', ms);
plot(freqdom(start_index:end), 100*DQE_RedCam.i1p717.a425(start_index:end), 'o-', 'Color', colors(2,:), 'LineWidth', lw, 'MarkerSize', ms);
plot(freqdom(start_index:end), 100*DQE_RedCam.i1p847.a425(start_index:end), 'd-', 'Color', colors(3,:), 'LineWidth', lw, 'MarkerSize', ms);
plot(freqdom(start_index:end), 100*DQE_RedCam.i1p923.a425(start_index:end), 'k-', 'Color', colors(4,:), 'LineWidth', lw, 'MarkerSize', ms);
plot(freqdom(start_index:end), 100*DQE_RedCam.i2p016.a425(start_index:end), 's-', 'Color', colors(5,:), 'LineWidth', lw, 'MarkerSize', ms);
legend(flip(red_int_gains), 'FontSize', 16);
legend boxoff
ylabel('DQE (%)');
xlabel('Frequency (mm^{-1})');
xlim([xmin,xmax]);
ylim([0,1.7]);
ax = gca;
ax.FontSize = 16; 
grid on;


subplot(1,2,2);
text(labelPos(1),labelPos(2),'b)', 'Units', 'normalized', 'Fontsize', 16)
hold on;
plot(freqdom(start_index:end), 100*DQE_LNCam.i1p848.a425(start_index:end), '^-', 'Color', colors(1,:), 'LineWidth', lw, 'MarkerSize', ms);
plot(freqdom(start_index:end), 100*DQE_LNCam.i1p922.a425(start_index:end), 'o-', 'Color', colors(2,:), 'LineWidth', lw, 'MarkerSize', ms);
plot(freqdom(start_index:end), 100*DQE_LNCam.i2p034.a425(start_index:end), 'd-', 'Color', colors(3,:), 'LineWidth', lw, 'MarkerSize', ms);
plot(freqdom(start_index:end), 100*DQE_LNCam.i2p135.a425(start_index:end), 'v-', 'Color', colors(4,:), 'LineWidth', lw, 'MarkerSize', ms);
plot(freqdom(start_index:end), 100*DQE_LNCam.i2p256.a425(start_index:end), 's-', 'Color', colors(5,:), 'LineWidth', lw, 'MarkerSize', ms);
legend(flip(LN_int_gains), 'FontSize', 16);
legend boxoff
ylabel('DQE (%)');
xlabel('Frequency (mm^{-1})');
xlim([xmin,xmax]);
ylim([0,1.7]);
ax = gca;
ax.FontSize = 16; 
grid on;

cd '/Users/danielalexander/Documents/Analysis/DQE/DQE Analysis/'
saveas(gcf, 'Figs/DQE_Int.png')

%% DQE - Compare Filtered/Unfiltered

figure('position', figSize2);
hold on;
plot(freqdom(start_index:end), 100*DQE_RedCam.i1p602.a425(start_index:end),       '^-',  'Color', colors(1,:), 'LineWidth', lw, 'MarkerSize', ms);
plot(freqdom(start_index:end), 100*DQE_RedCam_filt.i1p602.a425(start_index:end),  '^--', 'Color', colors(1,:), 'LineWidth', lw, 'MarkerSize', ms, 'MarkerFaceColor', colors(1,:));
plot(freqdom(start_index:end), 100*DQE_LNCam.i1p848.a425(start_index:end),        's-',  'Color', colors(2,:), 'LineWidth', lw, 'MarkerSize', ms);
plot(freqdom(start_index:end), 100*DQE_LNCam_filt.i1p848.a425(start_index:end),   's--', 'Color', colors(2,:), 'LineWidth', lw, 'MarkerSize', ms, 'MarkerFaceColor', colors(2,:));
ylabel('DQE (%)');
xlabel('Frequency (mm^{-1})');
xlim([xmin,xmax]);
ylim([0,1.7]);
legend('Gen3 Raw', 'Gen3 Filtered', 'Gen2+ Raw', 'Gen2+ Filtered');
legend boxoff
grid on;
ax = gca;
ax.FontSize = 16; 

cd '/Users/danielalexander/Documents/Analysis/DQE/DQE Analysis/'
saveas(gcf, 'Figs/DQE_Comp.png')

%% DF NPS

load('/Volumes/GoogleDrive/My Drive/Analysis/DQE/DQE Analysis/Results/DF_2DNPS_results.mat');

Fs = 1/pixelsize;
[Fx, Fy] = meshgrid((-Fs/2 + Fs/n):Fs/n:Fs/2, (-Fs/2 + Fs/n):Fs/n:Fs/2);


figure('position', figSize2);
% plot on large axes
surf(Fx, Fy, mean(NPS2D_RedCam_DF,3))% create smaller axes in top right, and plot on it
xlabel('f_x (mm^{-1})');
ylabel('f_y (mm^{-1})');
zlabel('NPS (mm^2)');
set(gca, 'FontSize', 16);
colormap(jet);
axes('Position',[.6 .6 .3 .3])
box on
imagesc(mean(NPS2D_RedCam_DF,3));
axis image; axis off;
colormap(jet);
cd '/Volumes/GoogleDrive/My Drive/Analysis/DQE/DQE Analysis/'
saveas(gcf, 'Figs/DF_NPS_Red.png')


figure('position', figSize2);
% plot on large axes
surf(Fx, Fy, mean(NPS2D_LNCam_DF,3))% create smaller axes in top right, and plot on it
xlabel('f_x (mm^{-1})');
ylabel('f_y (mm^{-1})');
zlabel('NPS (mm^2)');
set(gca, 'FontSize', 16);
colormap(jet);
axes('Position',[.6 .6 .3 .3])
box on
imagesc(mean(NPS2D_LNCam_DF,3));
axis image; axis off;
colormap(jet);
cd '/Volumes/GoogleDrive/My Drive/Analysis/DQE/DQE Analysis/'
saveas(gcf, 'Figs/DF_NPS_LN.png')


%% 2D NPS

load('/Volumes/GoogleDrive/My Drive/Analysis/DQE/DQE Analysis/Results/2DNPS_results.mat');

Fs = 1/pixelsize;
[Fx, Fy] = meshgrid((-Fs/2 + Fs/n):Fs/n:Fs/2, (-Fs/2 + Fs/n):Fs/n:Fs/2);

figure('position', figSize2);
% plot on large axes
surf(Fx, Fy, mean(NPS2D_RedCam.i1p602.a425,3))% create smaller axes in top right, and plot on it
xlabel('f_x (mm^{-1})');
ylabel('f_y (mm^{-1})');
zlabel('NPS (mm^2)');
set(gca, 'FontSize', 16);
colormap(jet);
axes('Position',[.6 .6 .3 .3])
box on
imagesc(mean(NPS2D_RedCam.i1p602.a425,3));
axis image; axis off;
colormap(jet);
cd '/Volumes/GoogleDrive/My Drive/Analysis/DQE/DQE Analysis/'
saveas(gcf, 'Figs/2D_NPS_Red.png')

figure('position', figSize2);
% plot on large axes
surf(Fx, Fy, mean(NPS2D_LNCam.i1p848.a425,3))% create smaller axes in top right, and plot on it
xlabel('f_x (mm^{-1})');
ylabel('f_y (mm^{-1})');
zlabel('NPS (mm^2)');
set(gca, 'FontSize', 16);
colormap(jet);
axes('Position',[.6 .6 .3 .3])
box on
imagesc(mean(NPS2D_LNCam.i1p848.a425,3));
axis image; axis off;
colormap(jet);
cd '/Volumes/GoogleDrive/My Drive/Analysis/DQE/DQE Analysis/'
saveas(gcf, 'Figs/2D_NPS_LN.png')


