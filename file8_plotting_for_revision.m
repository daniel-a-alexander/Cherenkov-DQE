%% Path

addpath(genpath('/Users/danielalexander/OneDrive - Dartmouth College/Matlab Functions/'));

%% Load

load('Results/revisionErrorAnalysis_LNCam.mat');
load('Results/revisionErrorAnalysis_RedCam.mat');
load('Results/revisionErrorAnalysis_filt.mat');

%% get standard devs, standard errors, avgs

n_frames = 159;

fields_ln = fieldnames(d_LNCam);
fields_rd = fieldnames(d_RedCam);

d_stats_LNCam = struct();
d_stats_RedCam = struct();

NPS1D_stats_LNCam = struct();
NPS1D_stats_RedCam = struct();


for i = 1:numel(fields_ln)
    sec_fields_ln = fieldnames(d_LNCam.(fields_ln{i}));
    sec_fields_rd = fieldnames(d_RedCam.(fields_rd{i}));
    for j = 1:numel(sec_fields_ln)
        % LN cam d
        d_stats_LNCam.(fields_ln{i}).(sec_fields_ln{j}).avg = mean(d_LNCam.(fields_ln{i}).(sec_fields_ln{j}),2);
        d_stats_LNCam.(fields_ln{i}).(sec_fields_ln{j}).stddev = std(d_LNCam.(fields_ln{i}).(sec_fields_ln{j}),0,2);
        d_stats_LNCam.(fields_ln{i}).(sec_fields_ln{j}).stderr = d_stats_LNCam.(fields_ln{i}).(sec_fields_ln{j}).stddev./sqrt(n_frames);
        
        % Red cam d
        d_stats_RedCam.(fields_rd{i}).(sec_fields_rd{j}).avg = mean(d_RedCam.(fields_rd{i}).(sec_fields_rd{j}),2);
        d_stats_RedCam.(fields_rd{i}).(sec_fields_rd{j}).stddev = std(d_RedCam.(fields_rd{i}).(sec_fields_rd{j}),0,2);
        d_stats_RedCam.(fields_rd{i}).(sec_fields_rd{j}).stderr = d_stats_RedCam.(fields_rd{i}).(sec_fields_rd{j}).stddev./sqrt(n_frames);
        
        % LN cam NPS
        NPS1D_stats_LNCam.(fields_ln{i}).(sec_fields_ln{j}).avg = mean(NPS1D_LNCam.(fields_ln{i}).(sec_fields_ln{j}),2);
        NPS1D_stats_LNCam.(fields_ln{i}).(sec_fields_ln{j}).stddev = std(NPS1D_LNCam.(fields_ln{i}).(sec_fields_ln{j}),0,2);
        NPS1D_stats_LNCam.(fields_ln{i}).(sec_fields_ln{j}).stderr = NPS1D_stats_LNCam.(fields_ln{i}).(sec_fields_ln{j}).stddev./sqrt(n_frames);
        
        % Red cam NPS
        NPS1D_stats_RedCam.(fields_rd{i}).(sec_fields_rd{j}).avg = mean(NPS1D_RedCam.(fields_rd{i}).(sec_fields_rd{j}),2);
        NPS1D_stats_RedCam.(fields_rd{i}).(sec_fields_rd{j}).stddev = std(NPS1D_RedCam.(fields_rd{i}).(sec_fields_rd{j}),0,2);
        NPS1D_stats_RedCam.(fields_rd{i}).(sec_fields_rd{j}).stderr = NPS1D_stats_RedCam.(fields_rd{i}).(sec_fields_rd{j}).stddev./sqrt(n_frames);
    end
end

%% Do it for filt data

NPS1D_stats_LNCam_filt = struct();

NPS1D_stats_LNCam_filt.i1p848.a425.avg = mean(NPS1D_LNCam_filt.i1p848.a425,2);
NPS1D_stats_LNCam_filt.i1p848.a425.stddev = std(NPS1D_LNCam_filt.i1p848.a425(5:end),0,2);
NPS1D_stats_LNCam_filt.i1p848.a425.stderr = NPS1D_stats_LNCam_filt.i1p848.a425.stddev./sqrt(n_frames);

NPS1D_stats_RedCam_filt = struct();

NPS1D_stats_RedCam_filt.i1p602.a425.avg = mean(NPS1D_RedCam_filt.i1p602.a425,2);
NPS1D_stats_RedCam_filt.i1p602.a425.stddev = std(NPS1D_RedCam_filt.i1p602.a425(5:end),0,2);
NPS1D_stats_RedCam_filt.i1p602.a425.stderr = NPS1D_stats_RedCam_filt.i1p602.a425.stddev./sqrt(n_frames);

%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% PLOTTING %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Initialize

load('/Users/danielalexander/OneDrive - Dartmouth College/Research/DQE/DQE Analysis/Results/DQE_all_results.mat', 'r');

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


%% NPS - Vary CMOS Gain

figure('position', figSize1);

subplot(1,2,1);
text(labelPos(1),labelPos(2),'a)', 'Units', 'normalized', 'Fontsize', 16)
hold on;
errorbar(freqdom(start_index:end), NPS1D_stats_RedCam.i1p602.a425.avg(start_index:end), ...
    NPS1D_stats_RedCam.i1p602.a425.stddev(start_index:end), '^-', 'Color', colors(1,:), 'LineWidth', lw, 'MarkerSize', ms);
errorbar(freqdom(start_index:end), NPS1D_stats_RedCam.i1p602.a400.avg(start_index:end), ...
    NPS1D_stats_RedCam.i1p602.a400.stddev(start_index:end), 'o-', 'Color', colors(2,:), 'LineWidth', lw, 'MarkerSize', ms);
errorbar(freqdom(start_index:end), NPS1D_stats_RedCam.i1p602.a375.avg(start_index:end), ...
    NPS1D_stats_RedCam.i1p602.a375.stddev(start_index:end), 'd-', 'Color', colors(3,:), 'LineWidth', lw, 'MarkerSize', ms);
errorbar(freqdom(start_index:end), NPS1D_stats_RedCam.i1p602.a350.avg(start_index:end), ...
    NPS1D_stats_RedCam.i1p602.a350.stddev(start_index:end), 'v-', 'Color', colors(4,:), 'LineWidth', lw, 'MarkerSize', ms);
errorbar(freqdom(start_index:end), NPS1D_stats_RedCam.i1p602.a325.avg(start_index:end), ...
    NPS1D_stats_RedCam.i1p602.a325.stddev(start_index:end), 's-', 'Color', colors(5,:), 'LineWidth', lw, 'MarkerSize', ms);
errorbar(freqdom(start_index:end), NPS1D_stats_RedCam.i1p602.a300.avg(start_index:end), ...
    NPS1D_stats_RedCam.i1p602.a300.stddev(start_index:end), '+-', 'Color', colors(6,:), 'LineWidth', lw, 'MarkerSize', ms);
errorbar(freqdom(start_index:end), NPS1D_stats_RedCam.i1p602.a250.avg(start_index:end), ...
    NPS1D_stats_RedCam.i1p602.a250.stddev(start_index:end), '*-', 'Color', colors(7,:), 'LineWidth', lw, 'MarkerSize', ms);
errorbar(freqdom(start_index:end), NPS1D_stats_RedCam.i1p602.a200.avg(start_index:end), ...
    NPS1D_stats_RedCam.i1p602.a200.stddev(start_index:end), 'p-', 'Color', colors(8,:), 'LineWidth', lw, 'MarkerSize', ms);
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
errorbar(freqdom(start_index:end), NPS1D_stats_LNCam.i1p848.a425.avg(start_index:end), ...
    NPS1D_stats_LNCam.i1p848.a425.stddev(start_index:end), '^-', 'Color', colors(1,:), 'LineWidth', lw, 'MarkerSize', ms);
errorbar(freqdom(start_index:end), NPS1D_stats_LNCam.i1p848.a400.avg(start_index:end), ...
    NPS1D_stats_LNCam.i1p848.a400.stddev(start_index:end), 'o-', 'Color', colors(2,:), 'LineWidth', lw, 'MarkerSize', ms);
errorbar(freqdom(start_index:end), NPS1D_stats_LNCam.i1p848.a375.avg(start_index:end), ...
    NPS1D_stats_LNCam.i1p848.a375.stddev(start_index:end), 'd-', 'Color', colors(3,:), 'LineWidth', lw, 'MarkerSize', ms);
errorbar(freqdom(start_index:end), NPS1D_stats_LNCam.i1p848.a350.avg(start_index:end), ...
    NPS1D_stats_LNCam.i1p848.a350.stddev(start_index:end), 'v-', 'Color', colors(4,:), 'LineWidth', lw, 'MarkerSize', ms);
errorbar(freqdom(start_index:end), NPS1D_stats_LNCam.i1p848.a325.avg(start_index:end), ...
    NPS1D_stats_LNCam.i1p848.a325.stddev(start_index:end), 's-', 'Color', colors(5,:), 'LineWidth', lw, 'MarkerSize', ms);
errorbar(freqdom(start_index:end), NPS1D_stats_LNCam.i1p848.a300.avg(start_index:end), ...
    NPS1D_stats_LNCam.i1p848.a300.stddev(start_index:end), '+-', 'Color', colors(6,:), 'LineWidth', lw, 'MarkerSize', ms);
errorbar(freqdom(start_index:end), NPS1D_stats_LNCam.i1p848.a250.avg(start_index:end), ...
    NPS1D_stats_LNCam.i1p848.a250.stddev(start_index:end), '*-', 'Color', colors(7,:), 'LineWidth', lw, 'MarkerSize', ms);
errorbar(freqdom(start_index:end), NPS1D_stats_LNCam.i1p848.a200.avg(start_index:end), ...
    NPS1D_stats_LNCam.i1p848.a200.stddev(start_index:end), 'p-', 'Color', colors(8,:), 'LineWidth', lw, 'MarkerSize', ms);
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

% cd '/Users/danielalexander/Documents/Analysis/DQE/DQE Analysis/'
saveas(gcf, 'Revised Figs/NPS_CMOS.png')

%% NPS - vary int gain

figure('position', figSize1);

subplot(1,2,1);
text(labelPos(1),labelPos(2),'a)', 'Units', 'normalized', 'Fontsize', 16)
hold on;
errorbar(freqdom(start_index:end), NPS1D_stats_RedCam.i1p602.a425.avg(start_index:end), ...
    NPS1D_stats_RedCam.i1p602.a425.stddev(start_index:end), '^-', 'Color', colors(1,:), 'LineWidth', lw, 'MarkerSize', ms);
errorbar(freqdom(start_index:end), NPS1D_stats_RedCam.i1p717.a425.avg(start_index:end), ...
    NPS1D_stats_RedCam.i1p717.a425.stddev(start_index:end), 'o-', 'Color', colors(2,:), 'LineWidth', lw, 'MarkerSize', ms);
errorbar(freqdom(start_index:end), NPS1D_stats_RedCam.i1p847.a425.avg(start_index:end), ...
    NPS1D_stats_RedCam.i1p847.a425.stddev(start_index:end), 'd-', 'Color', colors(3,:), 'LineWidth', lw, 'MarkerSize', ms);
errorbar(freqdom(start_index:end), NPS1D_stats_RedCam.i1p923.a425.avg(start_index:end), ...
    NPS1D_stats_RedCam.i1p923.a425.stddev(start_index:end), 'v-', 'Color', colors(4,:), 'LineWidth', lw, 'MarkerSize', ms);
errorbar(freqdom(start_index:end), NPS1D_stats_RedCam.i2p016.a425.avg(start_index:end), ...
    NPS1D_stats_RedCam.i2p016.a425.stddev(start_index:end), 's-', 'Color', colors(5,:), 'LineWidth', lw, 'MarkerSize', ms);
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
errorbar(freqdom(start_index:end), NPS1D_stats_LNCam.i1p848.a425.avg(start_index:end), ...
    NPS1D_stats_LNCam.i1p848.a425.stddev(start_index:end), '^-', 'Color', colors(1,:), 'LineWidth', lw, 'MarkerSize', ms);
errorbar(freqdom(start_index:end), NPS1D_stats_LNCam.i1p922.a425.avg(start_index:end), ...
    NPS1D_stats_LNCam.i1p922.a425.stddev(start_index:end), 'o-', 'Color', colors(2,:), 'LineWidth', lw, 'MarkerSize', ms);
errorbar(freqdom(start_index:end), NPS1D_stats_LNCam.i2p034.a425.avg(start_index:end), ...
    NPS1D_stats_LNCam.i2p034.a425.stddev(start_index:end), 'd-', 'Color', colors(3,:), 'LineWidth', lw, 'MarkerSize', ms);
errorbar(freqdom(start_index:end), NPS1D_stats_LNCam.i2p135.a425.avg(start_index:end), ...
    NPS1D_stats_LNCam.i2p135.a425.stddev(start_index:end), 'v-', 'Color', colors(4,:), 'LineWidth', lw, 'MarkerSize', ms);
errorbar(freqdom(start_index:end), NPS1D_stats_LNCam.i2p256.a425.avg(start_index:end), ...
    NPS1D_stats_LNCam.i2p256.a425.stddev(start_index:end), 's-', 'Color', colors(5,:), 'LineWidth', lw, 'MarkerSize', ms);
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

% cd '/Users/danielalexander/Documents/Analysis/DQE/DQE Analysis/'
saveas(gcf, 'Revised Figs/NPS_int.png')

%% NPS - Compare Filtered/Unfiltered

figure('position', figSize2);

hold on;
errorbar(freqdom(start_index:end), NPS1D_stats_RedCam.i1p602.a425.avg(start_index:end),...
    NPS1D_stats_RedCam.i1p602.a425.stddev(start_index:end), '^-',  'Color', colors(1,:), 'LineWidth', lw, 'MarkerSize', ms);
errorbar(freqdom(start_index:end), NPS1D_stats_RedCam_filt.i1p602.a425.avg(start_index:end),...
    NPS1D_stats_RedCam_filt.i1p602.a425.stddev(start_index:end), '^--', 'Color', colors(1,:), 'LineWidth', lw, 'MarkerSize', ms, 'MarkerFaceColor', colors(1,:));
errorbar(freqdom(start_index:end), NPS1D_stats_LNCam.i1p848.a425.avg(start_index:end),...
    NPS1D_stats_LNCam.i1p848.a425.stddev(start_index:end), 's-',  'Color', colors(2,:), 'LineWidth', lw, 'MarkerSize', ms);
errorbar(freqdom(start_index:end), NPS1D_stats_LNCam_filt.i1p848.a425.avg(start_index:end),...
     NPS1D_stats_LNCam_filt.i1p848.a425.stddev(start_index:end), 's--', 'Color', colors(2,:), 'LineWidth', lw, 'MarkerSize', ms, 'MarkerFaceColor', colors(2,:));
set(gca, 'YScale', 'log')
ylabel('NPS (mm^2)');
xlabel('Frequency (mm^{-1})');
xlim([xmin,xmax]);
legend('Gen3 Raw', 'Gen3 Filtered', 'Gen2+ Raw', 'Gen2+ Filtered');
legend boxoff
grid on;
ax = gca;
ax.FontSize = 16; 

% cd '/Users/danielalexander/Documents/Analysis/DQE/DQE Analysis/'
saveas(gcf, 'Revised Figs/NPS_Comp.png')

%% Load for DQE plotting

load('/Users/danielalexander/OneDrive - Dartmouth College/Research/DQE/DQE Analysis/Results/DQE_all_results.mat');

%% DQE - Vary CMOS Gain

figure('position', figSize1);

subplot(1,2,1);
text(labelPos(1),labelPos(2),'a)', 'Units', 'normalized', 'Fontsize', 16)
hold on;

errorbar(freqdom(start_index:end), 100*DQE_RedCam.i1p602.a425(start_index:end),...
    100*DQE_RedCam.i1p602.a425(start_index:end)./NPS1D_stats_RedCam.i1p602.a425.avg(start_index:end).*NPS1D_stats_RedCam.i1p602.a425.stddev(start_index:end),...
    '^-', 'Color', colors(1,:), 'LineWidth', lw, 'MarkerSize', ms);

errorbar(freqdom(start_index:end), 100*DQE_RedCam.i1p602.a400(start_index:end),...
    100*DQE_RedCam.i1p602.a400(start_index:end)./NPS1D_stats_RedCam.i1p602.a400.avg(start_index:end).*NPS1D_stats_RedCam.i1p602.a400.stddev(start_index:end),... 
    'o-', 'Color', colors(2,:), 'LineWidth', lw, 'MarkerSize', ms);

errorbar(freqdom(start_index:end), 100*DQE_RedCam.i1p602.a375(start_index:end),...
    100*DQE_RedCam.i1p602.a375(start_index:end)./NPS1D_stats_RedCam.i1p602.a375.avg(start_index:end).*NPS1D_stats_RedCam.i1p602.a375.stddev(start_index:end),... 
    'd-', 'Color', colors(3,:), 'LineWidth', lw, 'MarkerSize', ms);

errorbar(freqdom(start_index:end), 100*DQE_RedCam.i1p602.a350(start_index:end),...
    100*DQE_RedCam.i1p602.a350(start_index:end)./NPS1D_stats_RedCam.i1p602.a350.avg(start_index:end).*NPS1D_stats_RedCam.i1p602.a350.stddev(start_index:end),...
    'v-', 'Color', colors(4,:), 'LineWidth', lw, 'MarkerSize', ms);

errorbar(freqdom(start_index:end), 100*DQE_RedCam.i1p602.a325(start_index:end),...
    100*DQE_RedCam.i1p602.a325(start_index:end)./NPS1D_stats_RedCam.i1p602.a325.avg(start_index:end).*NPS1D_stats_RedCam.i1p602.a325.stddev(start_index:end),...
    's-', 'Color', colors(5,:), 'LineWidth', lw, 'MarkerSize', ms);

errorbar(freqdom(start_index:end), 100*DQE_RedCam.i1p602.a300(start_index:end),...
    100*DQE_RedCam.i1p602.a300(start_index:end)./NPS1D_stats_RedCam.i1p602.a300.avg(start_index:end).*NPS1D_stats_RedCam.i1p602.a300.stddev(start_index:end),...
    '+-', 'Color', colors(6,:), 'LineWidth', lw, 'MarkerSize', ms);

errorbar(freqdom(start_index:end), 100*DQE_RedCam.i1p602.a250(start_index:end),...
    100*DQE_RedCam.i1p602.a250(start_index:end)./NPS1D_stats_RedCam.i1p602.a250.avg(start_index:end).*NPS1D_stats_RedCam.i1p602.a250.stddev(start_index:end),...
    '*-', 'Color', colors(7,:), 'LineWidth', lw, 'MarkerSize', ms);

errorbar(freqdom(start_index:end), 100*DQE_RedCam.i1p602.a200(start_index:end),...
    100*DQE_RedCam.i1p602.a200(start_index:end)./NPS1D_stats_RedCam.i1p602.a200.avg(start_index:end).*NPS1D_stats_RedCam.i1p602.a200.stddev(start_index:end),...
    'p-', 'Color', colors(8,:), 'LineWidth', lw, 'MarkerSize', ms);

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

errorbar(freqdom(start_index:end), 100*DQE_LNCam.i1p848.a425(start_index:end),...
    100*DQE_LNCam.i1p848.a425(start_index:end)./NPS1D_stats_LNCam.i1p848.a425.avg(start_index:end).*NPS1D_stats_LNCam.i1p848.a425.stddev(start_index:end),...
    '^-', 'Color', colors(1,:), 'LineWidth', lw, 'MarkerSize', ms);

errorbar(freqdom(start_index:end), 100*DQE_LNCam.i1p848.a400(start_index:end),...
    100*DQE_LNCam.i1p848.a400(start_index:end)./NPS1D_stats_LNCam.i1p848.a400.avg(start_index:end).*NPS1D_stats_LNCam.i1p848.a400.stddev(start_index:end),...
    'o-', 'Color', colors(2,:), 'LineWidth', lw, 'MarkerSize', ms);

errorbar(freqdom(start_index:end), 100*DQE_LNCam.i1p848.a375(start_index:end),...
    100*DQE_LNCam.i1p848.a375(start_index:end)./NPS1D_stats_LNCam.i1p848.a375.avg(start_index:end).*NPS1D_stats_LNCam.i1p848.a375.stddev(start_index:end),...
    'd-', 'Color', colors(3,:), 'LineWidth', lw, 'MarkerSize', ms);

errorbar(freqdom(start_index:end), 100*DQE_LNCam.i1p848.a350(start_index:end),...
    100*DQE_LNCam.i1p848.a350(start_index:end)./NPS1D_stats_LNCam.i1p848.a350.avg(start_index:end).*NPS1D_stats_LNCam.i1p848.a350.stddev(start_index:end),...
    'v-', 'Color', colors(4,:), 'LineWidth', lw, 'MarkerSize', ms);

errorbar(freqdom(start_index:end), 100*DQE_LNCam.i1p848.a325(start_index:end),...
    100*DQE_LNCam.i1p848.a325(start_index:end)./NPS1D_stats_LNCam.i1p848.a325.avg(start_index:end).*NPS1D_stats_LNCam.i1p848.a325.stddev(start_index:end),...
    's-', 'Color', colors(5,:), 'LineWidth', lw, 'MarkerSize', ms);

errorbar(freqdom(start_index:end), 100*DQE_LNCam.i1p848.a300(start_index:end),...
    100*DQE_LNCam.i1p848.a300(start_index:end)./NPS1D_stats_LNCam.i1p848.a300.avg(start_index:end).*NPS1D_stats_LNCam.i1p848.a300.stddev(start_index:end),...
    '+-', 'Color', colors(6,:), 'LineWidth', lw, 'MarkerSize', ms);

errorbar(freqdom(start_index:end), 100*DQE_LNCam.i1p848.a250(start_index:end),...
    100*DQE_LNCam.i1p848.a250(start_index:end)./NPS1D_stats_LNCam.i1p848.a250.avg(start_index:end).*NPS1D_stats_LNCam.i1p848.a250.stddev(start_index:end),...
    '*-', 'Color', colors(7,:), 'LineWidth', lw, 'MarkerSize', ms);

errorbar(freqdom(start_index:end), 100*DQE_LNCam.i1p848.a200(start_index:end),...
    100*DQE_LNCam.i1p848.a200(start_index:end)./NPS1D_stats_LNCam.i1p848.a200.avg(start_index:end).*NPS1D_stats_LNCam.i1p848.a200.stddev(start_index:end),...
    'p-', 'Color', colors(8,:), 'LineWidth', lw, 'MarkerSize', ms);

legend(flip(cmos_gains), 'FontSize', 16, 'NumColumns',1);
legend boxoff
ylabel('DQE (%)');
xlabel('Frequency (mm^{-1})');
xlim([xmin,xmax]);
ylim([0,1.7]);
ax = gca;
ax.FontSize = 16; 
grid on;



saveas(gcf, 'Revised Figs/DQE_CMOS.png')


%% DQE - Vary Intensifier Gain

figure('position', figSize1);

subplot(1,2,1);
text(labelPos(1),labelPos(2),'a)', 'Units', 'normalized', 'Fontsize', 16)
hold on;

errorbar(freqdom(start_index:end), 100*DQE_RedCam.i1p602.a425(start_index:end),...
    100*DQE_RedCam.i1p602.a425(start_index:end)./NPS1D_stats_RedCam.i1p602.a425.avg(start_index:end).*NPS1D_stats_RedCam.i1p602.a425.stddev(start_index:end),...
    '^-', 'Color', colors(1,:), 'LineWidth', lw, 'MarkerSize', ms);

errorbar(freqdom(start_index:end), 100*DQE_RedCam.i1p717.a425(start_index:end),...
    100*DQE_RedCam.i1p717.a425(start_index:end)./NPS1D_stats_RedCam.i1p717.a425.avg(start_index:end).*NPS1D_stats_RedCam.i1p717.a425.stddev(start_index:end),...
    'o-', 'Color', colors(2,:), 'LineWidth', lw, 'MarkerSize', ms);

errorbar(freqdom(start_index:end), 100*DQE_RedCam.i1p847.a425(start_index:end),...
    100*DQE_RedCam.i1p847.a425(start_index:end)./NPS1D_stats_RedCam.i1p847.a425.avg(start_index:end).*NPS1D_stats_RedCam.i1p847.a425.stddev(start_index:end),...
    'd-', 'Color', colors(3,:), 'LineWidth', lw, 'MarkerSize', ms);

errorbar(freqdom(start_index:end), 100*DQE_RedCam.i1p923.a425(start_index:end),...
    100*DQE_RedCam.i1p923.a425(start_index:end)./NPS1D_stats_RedCam.i1p923.a425.avg(start_index:end).*NPS1D_stats_RedCam.i1p923.a425.stddev(start_index:end),...
    'k-', 'Color', colors(4,:), 'LineWidth', lw, 'MarkerSize', ms);

errorbar(freqdom(start_index:end), 100*DQE_RedCam.i2p016.a425(start_index:end),...
    100*DQE_RedCam.i2p016.a425(start_index:end)./NPS1D_stats_RedCam.i2p016.a425.avg(start_index:end).*NPS1D_stats_RedCam.i2p016.a425.stddev(start_index:end),...
    's-', 'Color', colors(5,:), 'LineWidth', lw, 'MarkerSize', ms);

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

errorbar(freqdom(start_index:end), 100*DQE_LNCam.i1p848.a425(start_index:end),...
    100*DQE_LNCam.i1p848.a425(start_index:end)./NPS1D_stats_LNCam.i1p848.a425.avg(start_index:end).*NPS1D_stats_LNCam.i1p848.a425.stddev(start_index:end),...
    '^-', 'Color', colors(1,:), 'LineWidth', lw, 'MarkerSize', ms);

errorbar(freqdom(start_index:end), 100*DQE_LNCam.i1p922.a425(start_index:end),...
    100*DQE_LNCam.i1p922.a425(start_index:end)./NPS1D_stats_LNCam.i1p922.a425.avg(start_index:end).*NPS1D_stats_LNCam.i1p922.a425.stddev(start_index:end),...
    'o-', 'Color', colors(2,:), 'LineWidth', lw, 'MarkerSize', ms);

errorbar(freqdom(start_index:end), 100*DQE_LNCam.i2p034.a425(start_index:end),...
    100*DQE_LNCam.i2p034.a425(start_index:end)./NPS1D_stats_LNCam.i2p034.a425.avg(start_index:end).*NPS1D_stats_LNCam.i2p034.a425.stddev(start_index:end),...
    'd-', 'Color', colors(3,:), 'LineWidth', lw, 'MarkerSize', ms);

errorbar(freqdom(start_index:end), 100*DQE_LNCam.i2p135.a425(start_index:end),...
    100*DQE_LNCam.i2p135.a425(start_index:end)./NPS1D_stats_LNCam.i2p135.a425.avg(start_index:end).*NPS1D_stats_LNCam.i2p135.a425.stddev(start_index:end),...
    'v-', 'Color', colors(4,:), 'LineWidth', lw, 'MarkerSize', ms);

errorbar(freqdom(start_index:end), 100*DQE_LNCam.i2p256.a425(start_index:end),...
    100*DQE_LNCam.i2p256.a425(start_index:end)./NPS1D_stats_LNCam.i2p256.a425.avg(start_index:end).*NPS1D_stats_LNCam.i2p256.a425.stddev(start_index:end),...
    's-', 'Color', colors(5,:), 'LineWidth', lw, 'MarkerSize', ms);

legend(flip(LN_int_gains), 'FontSize', 16);
legend boxoff
ylabel('DQE (%)');
xlabel('Frequency (mm^{-1})');
xlim([xmin,xmax]);
ylim([0,1.7]);
ax = gca;
ax.FontSize = 16; 
grid on;

saveas(gcf, 'Revised Figs/DQE_int.png')


%% DQE - Compare Filtered/Unfiltered

figure('position', figSize2);
hold on;

errorbar(freqdom(start_index:end), 100*DQE_RedCam.i1p602.a425(start_index:end),...
    100*DQE_RedCam.i1p602.a425(start_index:end)./NPS1D_stats_RedCam.i1p602.a425.avg(start_index:end).*NPS1D_stats_RedCam.i1p602.a425.stddev(start_index:end),...
    '^-',  'Color', colors(1,:), 'LineWidth', lw, 'MarkerSize', ms);

errorbar(freqdom(start_index:end), 100*DQE_RedCam_filt.i1p602.a425(start_index:end),...
    100*DQE_RedCam_filt.i1p602.a425(start_index:end)./NPS1D_stats_RedCam_filt.i1p602.a425.avg(start_index:end).*NPS1D_stats_RedCam_filt.i1p602.a425.stddev(start_index:end),...
    '^--', 'Color', colors(1,:), 'LineWidth', lw, 'MarkerSize', ms, 'MarkerFaceColor', colors(1,:));

errorbar(freqdom(start_index:end), 100*DQE_LNCam.i1p848.a425(start_index:end),...
    100*DQE_LNCam.i1p848.a425(start_index:end)./NPS1D_stats_LNCam.i1p848.a425.avg(start_index:end).*NPS1D_stats_LNCam.i1p848.a425.stddev(start_index:end),...
    's-',  'Color', colors(2,:), 'LineWidth', lw, 'MarkerSize', ms);

errorbar(freqdom(start_index:end), 100*DQE_LNCam_filt.i1p848.a425(start_index:end),...
    100*DQE_LNCam_filt.i1p848.a425(start_index:end)./NPS1D_stats_LNCam_filt.i1p848.a425.avg(start_index:end).*NPS1D_stats_LNCam_filt.i1p848.a425.stddev(start_index:end),...
    's--', 'Color', colors(2,:), 'LineWidth', lw, 'MarkerSize', ms, 'MarkerFaceColor', colors(2,:));

ylabel('DQE (%)');
xlabel('Frequency (mm^{-1})');
xlim([xmin,xmax]);
ylim([0,1.7]);
legend('Gen3 Raw', 'Gen3 Filtered', 'Gen2+ Raw', 'Gen2+ Filtered');
legend boxoff
grid on;
ax = gca;
ax.FontSize = 16; 

saveas(gcf, 'Revised Figs/DQE_Comp.png')
