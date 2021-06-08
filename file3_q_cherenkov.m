cd '/Users/danielalexander/Documents/Dartmouth/DQE/DQE Analysis'
addpath(genpath('/Users/danielalexander/Documents/Matlab Functions'));

%%

cts_per_ph = 36;

ch_avg_frame = double(mean(read_dovi('2019-09-30 20-47-54-907/meas_s1_cam0.dovi'), 3));

%%
% 
figure;
imagesc(ch_avg_frame); colorbar; colormap(jet); axis off;
% rect = getrect;

rect = [633, 583, 317, 245];
roi = imcrop(ch_avg_frame, rect);

%% Get pixel size

%  imtool(mean(RedCam.i1p602.a425, 3),[]);
l = 80;
pixelsize = 30/l;

%%

mean_cts = mean(roi(:));

mean_ph = mean_cts / cts_per_ph;

mean_ph_per_mm = mean_ph / (pixelsize^2);

q_cherenkov = mean_ph_per_mm;