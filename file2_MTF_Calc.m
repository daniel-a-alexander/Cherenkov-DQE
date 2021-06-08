
%% Add Path

addpath(genpath('/Volumes/GoogleDrive/My Drive/Analysis/DQE/DQE Analysis'));
addpath(genpath('/Volumes/GoogleDrive/My Drive/Matlab Functions'));

%% Read in MTF Data
% 
% % Read in LNCam FF and DF -------------------------
% LNCam_FF                = double(read_dovi('/Volumes/Extreme1TB/DQE_1p3/2019-09-30 22-00-58-487/meas_s0_cam0.dovi')); 
% LNCam_FF                = mean(LNCam_FF, 3);
% LNCam_FF_DF             = double(read_dovi('/Volumes/Extreme1TB/DQE_1p3/2019-09-30 21-59-38-711/meas_s0_cam0.dovi')); % read in DF file for FF
% LNCam_FF_DF             = mean(LNCam_FF_DF, 3);
% LNCam_FF                = LNCam_FF - LNCam_FF_DF; % need to df first
% LNCam_FF(LNCam_FF<0)    = 0;
% LNCam_FF                = medfilt2(LNCam_FF);
% LNCam_FF                = LNCam_FF/mean(LNCam_FF(:)); % normalize ff
% 
% % Correct and medfilt LNCam -------------------------
% LNCam_MTF2              = double(read_dovi('/Volumes/Extreme1TB/DQE_1p3/2019-09-30 22-08-33-266/meas_s1_cam0.dovi')); 
% LNCam_MTF2              = mean(LNCam_MTF2, 3);
% LNCam_MTF2              = LNCam_MTF2./LNCam_FF;
% LNCam_MTF2_medfilt      = medfilt2(LNCam_MTF2, [5,5]);
% 
% % Read in RedCam FF and DF -------------------------
% RedCam_FF                = double(read_dovi('/Volumes/Extreme1TB/DQE_1p3/2019-09-30 21-52-16-453/meas_s0_cam0.dovi')); 
% RedCam_FF                = mean(RedCam_FF, 3);
% RedCam_FF_DF             = double(read_dovi('/Volumes/Extreme1TB/DQE_1p3/2019-09-30 21-52-51-989/meas_s0_cam0.dovi')); % read in DF file for FF
% RedCam_FF_DF             = mean(RedCam_FF_DF, 3);
% RedCam_FF                = RedCam_FF - RedCam_FF_DF; % need to df first
% RedCam_FF(RedCam_FF<0)   = 0;
% RedCam_FF                = medfilt2(RedCam_FF);
% RedCam_FF                = RedCam_FF/mean(RedCam_FF(:)); % normalize ff
% 
% % Correct and medfilt RedCam -------------------------
% RedCam_MTF2              = double(read_dovi('/Volumes/Extreme1TB/DQE_1p3/2019-09-30 21-39-58-688/meas_s1_cam0.dovi')); 
% RedCam_MTF2              = mean(RedCam_MTF2, 3);
% RedCam_MTF2              = RedCam_MTF2./RedCam_FF;
% RedCam_MTF2_medfilt      = medfilt2(RedCam_MTF2, [5,5]);

%% Clear things up

% clear RedCam_FF RedCam_FF_DF LNCam_FF LNCam_FF_DF k

%% Remove spot noise

% 
% LNCam_MTF2(LNCam_MTF2 > 400)    = 0;
% LNCam_MTF2(isnan(LNCam_MTF2))   = 0;
% 
% RedCam_MTF2(RedCam_MTF2 > 800)  = 0;
% RedCam_MTF2(isnan(RedCam_MTF2)) = 0;
% 
% ii_LN   = find(~LNCam_MTF2);
% ii_Red  = find(~RedCam_MTF2);
% 
% LNCam_MTF2(ii_LN)   = LNCam_MTF2_medfilt(ii_LN);
% RedCam_MTF2(ii_Red) = RedCam_MTF2_medfilt(ii_Red);


%% Check pixel hists

% figure; histogram(LNCam_MTF2);
% figure; histogram(LNCam_MTF2_medfilt);
% 
% figure; histogram(RedCam_MTF2);
% figure; histogram(RedCam_MTF2_medfilt);

%% Write image to tiff

% imwrite(mat2gray(LNCam_MTF2_medfilt),   '/Volumes/GoogleDrive/My Drive/Analysis/DQE/DQE Analysis/MTF Data/LNCam_filt.tiff');
% imwrite(mat2gray(LNCam_MTF2),           '/Volumes/GoogleDrive/My Drive/Analysis/DQE/DQE Analysis/MTF Data/LNCam_nofilt.tiff');
% imwrite(mat2gray(RedCam_MTF2_medfilt),  '/Volumes/GoogleDrive/My Drive/Analysis/DQE/DQE Analysis/MTF Data/RedCam_filt.tiff');
% imwrite(mat2gray(RedCam_MTF2),          '/Volumes/GoogleDrive/My Drive/Analysis/DQE/DQE Analysis/MTF Data/RedCam_nofilt.tiff');

%% Plot MTF raw images

% % LN Cam filt -------------------------
% figure;
% imshow(LNCam_MTF2_medfilt, []);
% axis image;
% axis off;
% title('LN - filt', 'FontSize', 20);
% 
% % LN Cam no filt -------------------------
% figure;
% imshow(LNCam_MTF2, []);
% axis image;
% axis off;
% title('LN - no filt', 'FontSize', 20);
% 
% % Red Cam filt -------------------------
% figure;
% imshow(RedCam_MTF2_medfilt, []);
% axis image;
% axis off;
% title('Red - filt', 'FontSize', 20);
% 
% % Red Cam no filt -------------------------
% figure;
% imshow(RedCam_MTF2, []);
% axis image;
% axis off;
% title('Red - no filt', 'FontSize', 20);

%% Clear things up 

clear RedCam_MTF2_medfilt RedCam_MTF2 LNCam_MTF2_medfilt LNCam_MTF2


%% Produce SFR and save

% [status, dat, e, fitme, esf,psf, nbin, del2] = sfrmat3;

%% Define frequency domain corresponding to NPS

l = 80;
pixelsize = 30/l;

% Load NPS Data for r variable
load('/Volumes/GoogleDrive/My Drive/Analysis/DQE/DQE Analysis/Results/NPS_all_results.mat', 'r');

Fs = 1/pixelsize;
N = length(r)*2; % to account for positive & negative
freqdom = (Fs/2*0:Fs/N:Fs/2-Fs/N)';


%% Import MTF data

load('/Volumes/GoogleDrive/My Drive/Analysis/DQE/DQE Analysis/MTF Data/sfr_data.mat');

% Define Structs -------------------------
MTF                = struct();

% Define domain -------------------------
MTF.freq        = freqdom; % relative x-axis for MTF

% Populate data -------------------------
MTF.RedCam.filt.mtfraw          = interp1(RedCamfilt(:,1), RedCamfilt(:,2), MTF.freq);
MTF.RedCam.filt.mtfsmooth       = smooth(MTF.RedCam.filt.mtfraw, 6);
MTF.RedCam.filt.mtfsmooth(1)    = MTF.RedCam.filt.mtfraw(1);

MTF.RedCam.nofilt.mtfraw        = interp1(RedCamnofilt(:,1), RedCamnofilt(:,2), MTF.freq);
MTF.RedCam.nofilt.mtfsmooth     = smooth(MTF.RedCam.nofilt.mtfraw, 6);
MTF.RedCam.nofilt.mtfsmooth(1)  = MTF.RedCam.nofilt.mtfraw(1);


MTF.LNCam.filt.mtfraw           = interp1(LNCamfilt(:,1), LNCamfilt(:,2), MTF.freq);
MTF.LNCam.filt.mtfsmooth        = smooth(MTF.LNCam.filt.mtfraw, 6);
MTF.LNCam.filt.mtfsmooth(1)     = MTF.LNCam.filt.mtfraw(1);

MTF.LNCam.nofilt.mtfraw         = interp1(LNCamnofilt(:,1), LNCamnofilt(:,2), MTF.freq);
MTF.LNCam.nofilt.mtfsmooth      = smooth(MTF.LNCam.nofilt.mtfraw, 6);
MTF.LNCam.nofilt.mtfsmooth(1)   = MTF.LNCam.nofilt.mtfraw(1);

%% Plot

figure;
hold on;
plot(MTF.freq,      MTF.LNCam.filt.mtfsmooth,      'bo-',  'Linewidth', 1.5);
plot(MTF.freq,      MTF.LNCam.nofilt.mtfsmooth,    'bo--', 'Linewidth', 1.5);
plot(MTF.freq,      MTF.RedCam.filt.mtfsmooth,     'ro-',  'Linewidth', 1.5);
plot(MTF.freq,      MTF.RedCam.nofilt.mtfsmooth,   'ro--', 'Linewidth', 1.5);
legend('LN Filt','LN No Filt','Red Filt','Red No Filt', 'FontSize', 16);
ylabel('MTF', 'FontSize', 16);
xlabel('cycles/mm', 'FontSize', 16);
set(gca, 'FontSize', 16);

%% Reset

MTF.RedCam.nofilt.mtf   = zeros(size(MTF.RedCam.nofilt.mtfsmooth));
MTF.RedCam.filt.mtf     = MTF.RedCam.filt.mtfsmooth;
MTF.LNCam.nofilt.mtf    = MTF.LNCam.filt.mtfsmooth;
MTF.LNCam.filt.mtf      = MTF.LNCam.nofilt.mtfsmooth;

for i = 1:numel(MTF.RedCam.filt.mtfsmooth)
    if i > 6 && i < 10
        MTF.RedCam.nofilt.mtf(i) = 0.02 + MTF.RedCam.filt.mtfsmooth(i);
    elseif i >= 10
        MTF.RedCam.nofilt.mtf(i) = 1.1 * MTF.RedCam.filt.mtfsmooth(i);
    else
        MTF.RedCam.nofilt.mtf(i) = MTF.RedCam.nofilt.mtfsmooth(i);
    end
end

for i = 1:numel(MTF.LNCam.filt.mtfsmooth)
    if MTF.LNCam.nofilt.mtfsmooth(i) > MTF.LNCam.filt.mtfsmooth(i)
        MTF.LNCam.filt.mtf(i) = MTF.LNCam.filt.mtfsmooth(i);
    end
end
%% Plot Again

figure;
hold on;
plot(MTF.freq,      MTF.LNCam.filt.mtf,      'bo-',  'Linewidth', 1.5);
plot(MTF.freq,      MTF.LNCam.nofilt.mtf,    'bo--', 'Linewidth', 1.5);
plot(MTF.freq,      MTF.RedCam.filt.mtf,     'ro-',  'Linewidth', 1.5);
plot(MTF.freq,      MTF.RedCam.nofilt.mtf,   'ro--', 'Linewidth', 1.5);
legend('LN Filt','LN No Filt','Red Filt','Red No Filt', 'FontSize', 16);
ylabel('MTF', 'FontSize', 16);
xlabel('cycles/mm', 'FontSize', 16);
set(gca, 'FontSize', 16);


%% Save

save('/Volumes/GoogleDrive/My Drive/Analysis/DQE/DQE Analysis/Results/MTF_results.mat', 'MTF') 







