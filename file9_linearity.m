%% Add Path

% addpath(genpath('C:\Users\danal\OneDrive - Dartmouth College\Research\DQE\DQE Analysis'));
% addpath(genpath('C:\Users\danal\OneDrive - Dartmouth College\Matlab Functions'));
% 
% 
% cd 'C:\Users\danal\OneDrive - Dartmouth College\Research\DQE\DQE Analysis'

addpath(genpath('/Users/danielalexander/OneDrive - Dartmouth College/Research/DQE/DQE Analysis'));
addpath(genpath('/Users/danielalexander/OneDrive - Dartmouth College/Matlab Functions'));


cd '/Users/danielalexander/OneDrive - Dartmouth College/Research/DQE/DQE Analysis'


%% Read

% datapath = 'C:\Users\danal\OneDrive - Dartmouth College\Research\DQE\DQE Analysis\Linearity Data\data_071320';
datapath = '/Users/danielalexander/OneDrive - Dartmouth College/Research/DQE/DQE Analysis/Linearity Data/data_071320'

data_raw = struct();

data_raw.redcam.imbkg = read_dovi(fullfile(datapath,'2020-07-13 19-30-13-807/meas_s0_cam0.dovi'));
data_raw.redcam.im01 = read_dovi(fullfile(datapath,'2020-07-13 19-31-46-378/meas_s0_cam0.dovi'));
data_raw.redcam.im02 = read_dovi(fullfile(datapath,'2020-07-13 19-32-13-480/meas_s0_cam0.dovi'));
data_raw.redcam.im03 = read_dovi(fullfile(datapath,'2020-07-13 19-32-52-376/meas_s0_cam0.dovi'));
data_raw.redcam.im04 = read_dovi(fullfile(datapath,'2020-07-13 19-33-16-215/meas_s0_cam0.dovi'));
data_raw.redcam.im05 = read_dovi(fullfile(datapath,'2020-07-13 19-33-50-195/meas_s0_cam0.dovi'));
data_raw.redcam.im06 = read_dovi(fullfile(datapath,'2020-07-13 19-34-14-932/meas_s0_cam0.dovi'));
data_raw.redcam.im07 = read_dovi(fullfile(datapath,'2020-07-13 19-34-38-577/meas_s0_cam0.dovi'));
data_raw.redcam.im08 = read_dovi(fullfile(datapath,'2020-07-13 19-35-29-341/meas_s0_cam0.dovi'));
data_raw.redcam.im09 = read_dovi(fullfile(datapath,'2020-07-13 19-35-53-691/meas_s0_cam0.dovi'));

data_raw.lncam.imbkg = read_dovi(fullfile(datapath,'2020-07-13 19-45-13-143/meas_s0_cam0.dovi'));
data_raw.lncam.im01 = read_dovi(fullfile(datapath,'2020-07-13 19-46-01-971/meas_s0_cam0.dovi'));
data_raw.lncam.im02 = read_dovi(fullfile(datapath,'2020-07-13 19-46-45-347/meas_s0_cam0.dovi'));
data_raw.lncam.im03 = read_dovi(fullfile(datapath,'2020-07-13 19-47-45-489/meas_s0_cam0.dovi'));
data_raw.lncam.im04 = read_dovi(fullfile(datapath,'2020-07-13 19-48-08-900/meas_s0_cam0.dovi'));
data_raw.lncam.im05 = read_dovi(fullfile(datapath,'2020-07-13 19-49-05-810/meas_s0_cam0.dovi'));
data_raw.lncam.im06 = read_dovi(fullfile(datapath,'2020-07-13 19-49-44-505/meas_s0_cam0.dovi'));
data_raw.lncam.im07 = read_dovi(fullfile(datapath,'2020-07-13 19-50-09-050/meas_s0_cam0.dovi'));
data_raw.lncam.im08 = read_dovi(fullfile(datapath,'2020-07-13 19-50-35-700/meas_s0_cam0.dovi'));
data_raw.lncam.im09 = read_dovi(fullfile(datapath,'2020-07-13 19-50-56-242/meas_s0_cam0.dovi'));

%% bkg subtract

fields = fieldnames(data_raw.lncam);

redcambkg_avg = mean(data_raw.redcam.imbkg,3,'native'); % this takes a bit of time for native output
lncambkg_avg = mean(data_raw.lncam.imbkg,3,'native');

data = struct();

for i=2:numel(fields)
    for j=1:50
        data.redcam.(fields{i})(:,:,j) = data_raw.redcam.(fields{i})(:,:,j) - redcambkg_avg;
        data.lncam.(fields{i})(:,:,j) = data_raw.lncam.(fields{i})(:,:,j) - lncambkg_avg;
        disp(['frame', num2str(j)])
    end
end
    
%%  Analyze

imagesc(sum(data_raw.redcam.im09, 3));
redcam_rect = getrect();
    
imagesc(sum(data_raw.lncam.im09, 3));
lncam_rect = getrect();
    
data_cropped = struct();
data_avg = struct();
    
counts_lncam    = zeros(9,1);
counts_redcam   = zeros(9,1);
std_lncam       = zeros(9,1);
std_redcam      = zeros(9,1);

fields = fieldnames(data.lncam);
for i=1:numel(fields)
    for j=1:50
        data_cropped.redcam.(fields{i})(:,:,j)  = imcrop(data.redcam.(fields{i})(:,:,j),    redcam_rect);
        temp = data_cropped.redcam.(fields{i})(:,:,j);
        data_avg.redcam.(fields{i})(j) = mean(temp(:));
        
        data_cropped.lncam.(fields{i})(:,:,j)   = imcrop(data.lncam.(fields{i})(:,:,j),     lncam_rect);
        temp = data_cropped.lncam.(fields{i})(:,:,j);
        data_avg.lncam.(fields{i})(j) = mean(temp(:));
        
        disp(['frame', num2str(j)])
    end
    counts_redcam(i)    = mean(data_avg.redcam.(fields{i})(:));
    counts_lncam(i)     = mean(data_avg.lncam.(fields{i})(:));
    std_redcam(i)       = std(data_avg.redcam.(fields{i})(:));
    std_lncam(i)        = std(data_avg.lncam.(fields{i})(:));
end
    
%% Plotting

load('/Users/danielalexander/OneDrive - Dartmouth College/Research/DQE/DQE Analysis/Results/linearity_results.mat');
colors = linspecer(2);

power_redcam = [16.9, 33.2, 58.5, 102.5, 221, 441.5, 639.5, 836, 940];
power_lncam = [16.85, 33.2, 58.7, 102.4, 223, 442, 638, 835, 941];

xax = 1:1000;

figSize2 = [500,500,500,420];
ms = 6; % markersize
lw = 1; % Linewidth


figure('position', figSize2);

hold on
errorbar(power_redcam, counts_redcam, std_redcam, 'x', 'Color', colors(2,:), 'LineWidth', lw, 'MarkerSize', ms);
plot(xax, feval(fittedmodel_redcam, xax), '-', 'Color', colors(2,:), 'LineWidth', lw);
errorbar(power_lncam, counts_lncam, std_lncam, '*', 'Color', colors(1,:), 'LineWidth', lw, 'MarkerSize', ms);
plot(xax, feval(fittedmodel_lncam, xax), '-', 'Color', colors(1,:), 'LineWidth', lw);


ax = gca;
ax.FontSize = 16; 
% set(gca, 'YScale', 'log')
% set(gca, 'XScale', 'log')
legend('Gen3 data', 'Gen3 fit', 'Gen2+', 'Gen2+ fit', 'Location', 'Northwest');
% legend boxoff
xlabel('Power (µW)')
ylabel('Counts')
grid on
    
    
    
    
    
    
    
    

