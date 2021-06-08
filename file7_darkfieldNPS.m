
%% Add Path

addpath(genpath('/Volumes/GoogleDrive/My Drive/Analysis/DQE/DQE Analysis'));
addpath(genpath('/Volumes/GoogleDrive/My Drive/Matlab Functions'));

%% Read Data

datapath = '/Volumes/Extreme1TB/DQE_1p3/';

% Read in LNCam and DF ----------------------------------------------------------------------------------------------------
LNCam_DF             = read_dovi(fullfile(datapath, '2019-09-30 21-59-38-711/meas_s0_cam0.dovi')); 

% Read in RedCam and DF ----------------------------------------------------------------------------------------------------
RedCam_DF             = read_dovi(fullfile(datapath, '2019-09-30 21-52-51-989/meas_s0_cam0.dovi'));

%% Get pixel size

% checkerboard size is 30 mm

%  imtool(mean(LNCam_filt.i1p848.a425, 3),[]);
l = 80;
pixelsize = 30/l;

%% Get cropping rectangle from mean image if needed

LNCam_rect = [600, 500, 399, 399];
RedCam_rect = [600, 500, 399, 399];

%% Crop Images

frames = size(RedCam_DF,3);

A_Red = zeros(400, 400, frames);
A_LN = zeros(400, 400, frames);

for k=1:frames
    A_Red(:,:,k) = imcrop(RedCam_DF(:,:,k), RedCam_rect); % crop and assign
    A_LN(:,:,k) = imcrop(LNCam_DF(:,:,k), LNCam_rect); % crop and assign
end

A_RedCam_DF = A_Red;  % reassign new array to struct
A_LNCam_DF = A_LN;  % reassign new array to struct

%% NPS Calc Red Cam


m = 400; % size of image
n = 100; % size of each roi

NPS2D_RedCam_DF = zeros(n,n, frames); % preallocate size for matrix (contains 2D array for each frame)

for k=1:frames

    A_Red   = A_RedCam_DF(:,:,k); % currnt frame for nps
    
    FFTstack = zeros(n,n,length(1:n/2:m-n/2)^2); % this is the FFT for each ROI

    z = 1; % ROI index

    for g = 1:n/2:m-n/2
        for h = 1:n/2:m-n/2
            rect_gh = [g, h, n-1, n-1]; 
            Asub = imcrop(A_Red, rect_gh); % get ROI from cropped image 
            [X,Y]=meshgrid(1:n, 1:n); % define meshgrid to fit over
            [xData, yData, zData] = prepareSurfaceData( X, Y, Asub );
            ft = fittype( 'poly22' ); % Set up fittype and options
            [fitresult, gof] = fit( [xData, yData], zData, ft ); % Fit model to data
            S = fitresult(X,Y); % here's our fit result
            AsubC = (Asub - S); % we use this fit to correct for trends
            FFT = fft2(AsubC); % calc 2D FFT
            FFTstack(:,:,z) = abs(fftshift(FFT)).^2; % shift, scale amplitude, and calc mod squared and save to array
            z = z+1;
        end
    end

    %
    NPS = (pixelsize)^2*(1/n^2)*mean(FFTstack, 3);
    
    NPS2D_RedCam_DF(:,:,k) = NPS;
    
    disp([num2str(100*k/frames), ' % done']);

end


%% NPS Calc LN Cam

m = 400; % size of image
n = 100; % size of each roi

NPS2D_LNCam_DF = zeros(n,n, frames); % preallocate size for matrix (contains 2D array for each frame)

for k=1:frames
    
    A_LN    = A_LNCam_DF(:,:,k);
    
    FFTstack = zeros(n,n,length(1:n/2:m-n/2)^2); % this is the FFT for each ROI

    z = 1; % ROI index

    for g = 1:n/2:m-n/2
        for h = 1:n/2:m-n/2
            rect_gh = [g, h, n-1, n-1]; 
            Asub = imcrop(A_LN, rect_gh); % get ROI from cropped image 
            [X,Y]=meshgrid(1:n, 1:n); % define meshgrid to fit over
            [xData, yData, zData] = prepareSurfaceData( X, Y, Asub );
            ft = fittype( 'poly22' ); % Set up fittype and options
            [fitresult, gof] = fit( [xData, yData], zData, ft ); % Fit model to data
            S = fitresult(X,Y); % here's our fit result
            AsubC = (Asub - S); % we use this fit to correct for trends
            FFT = fft2(AsubC); % calc 2D FFT
            FFTstack(:,:,z) = abs(fftshift(FFT)).^2; % shift, scale amplitude, and calc mod squared and save to array
            z = z+1;
        end
    end

    %
    NPS = (pixelsize)^2*(1/n^2)*mean(FFTstack, 3);

    NPS2D_LNCam_DF(:,:,k) = NPS; % assign array to stuct of means
    
    disp([num2str(100*k/frames), ' % done']);

end


%% Clear

clearvars -except NPS2D_LNCam_DF NPS2D_RedCam_DF pixelsize n

%% Plot NPS


Fs = 1/pixelsize;
[Fx, Fy] = meshgrid((-Fs/2 + Fs/n):Fs/n:Fs/2, (-Fs/2 + Fs/n):Fs/n:Fs/2);


figure;
title('Gen3');
% plot on large axes
surf(Fx, Fy, mean(NPS2D_RedCam_DF,3))% create smaller axes in top right, and plot on it
xlabel('F_x (mm^{-1})');
ylabel('F_y (mm^{-1})');
zlabel('NPS (mm^2)');
set(gca, 'FontSize', 16);
% set(gca, 'ZScale', 'log')
% set(gca,'xtick',[],'ytick',[], 'ztick', []);
colormap(jet);
axes('Position',[.6 .6 .3 .3])
box on
imagesc(mean(NPS2D_RedCam_DF,3));
% xticks([0 25 50 75 100])
% xticklabels({num2str(-Fs/4),num2str(-Fs/4),'0',num2str(Fs/4),num2str(Fs/2)})
% yticks([0 25 50 75])
% yticklabels({num2str(Fs/4),num2str(Fs/4),'0',num2str(-Fs/4)})
axis image; axis off;
colormap(jet);

figure;
title('Gen2+');
% plot on large axes
surf(Fx, Fy, mean(NPS2D_LNCam_DF,3))% create smaller axes in top right, and plot on it
xlabel('F_x (mm^{-1})');
ylabel('F_y (mm^{-1})');
zlabel('NPS (mm^2)');
set(gca, 'FontSize', 16);
% set(gca, 'ZScale', 'log')
% set(gca,'xtick',[],'ytick',[], 'ztick', []);
colormap(jet);
axes('Position',[.6 .6 .3 .3])
box on
imagesc(mean(NPS2D_LNCam_DF,3));
% xticks([0 25 50 75 100])
% xticklabels({num2str(-Fs/4),num2str(-Fs/4),'0',num2str(Fs/4),num2str(Fs/2)})
% yticks([0 25 50 75])
% yticklabels({num2str(Fs/4),num2str(Fs/4),'0',num2str(-Fs/4)})
axis image; axis off;
colormap(jet);


%%








