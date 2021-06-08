%% Path

addpath(genpath('/Users/danielalexander/OneDrive - Dartmouth College/Matlab Functions/'));
datapath = '/Volumes/Extreme1TB/DQE_1p3/';

%% Data Read

% Read in LNCam FF and DF ----------------------------------------------------------------------------------------------------
LNCam_FF                = read_dovi(fullfile(datapath, '2019-09-30 22-00-58-487/meas_s0_cam0.dovi')); 
LNCam_FF                = mean(LNCam_FF, 3);
LNCam_FF_DF             = read_dovi(fullfile(datapath, '2019-09-30 21-59-38-711/meas_s0_cam0.dovi')); % read in DF file for FF
LNCam_FF_DF             = mean(LNCam_FF_DF, 3);
LNCam_FF                = LNCam_FF - LNCam_FF_DF; % need to df first
LNCam_FF(LNCam_FF<0)    = 0;
LNCam_FF                = medfilt2(LNCam_FF);
LNCam_FF                = LNCam_FF/mean(LNCam_FF(:)); % normalize ff

% Read in RedCam FF and DF ----------------------------------------------------------------------------------------------------
RedCam_FF                = read_dovi(fullfile(datapath, '2019-09-30 21-52-16-453/meas_s0_cam0.dovi')); 
RedCam_FF                = mean(RedCam_FF, 3);
RedCam_FF_DF             = read_dovi(fullfile(datapath, '2019-09-30 21-52-51-989/meas_s0_cam0.dovi')); % read in DF file for FF
RedCam_FF_DF             = mean(RedCam_FF_DF, 3);
RedCam_FF                = RedCam_FF - RedCam_FF_DF; % need to df first
RedCam_FF(RedCam_FF<0)   = 0;
RedCam_FF                = medfilt2(RedCam_FF);
RedCam_FF                = RedCam_FF/mean(RedCam_FF(:)); % normalize ff

% Define structs ----------------------------------------------------------------------------------------------------
RedCam_filt = struct;
LNCam_filt = struct;

% Read filtered stacks ----------------------------------------------------------------------------------------------------
RedCam_filt.i1p602.a425  = double(read_dovi('/Users/danielalexander/OneDrive - Dartmouth College/Research/DQE/DQE Analysis/temp_cdose_db_filt/data/DQE_1p3/2019-09-30 20-48-49-269/meas_s1_cam0.dovi'));
LNCam_filt.i1p848.a425   = double(read_dovi('/Users/danielalexander/OneDrive - Dartmouth College/Research/DQE/DQE Analysis/temp_cdose_db_filt/data/DQE_1p3/2019-09-30 19-53-44-158/meas_s1_cam0.dovi'));

% FF correct ----------------------------------------------------------------------------------------------------
% RedCam_filt.i1p602.a425     = double(RedCam_filt.i1p602.a425)./RedCam_FF;
% LNCam_filt.i1p848.a425      = double(LNCam_filt.i1p848.a425)./LNCam_FF;

frames = size(RedCam_filt.i1p602.a425,3);
for k=1:frames
    RedCam_filt.i1p602.a425(:,:,k) = RedCam_filt.i1p602.a425(:,:,k)./RedCam_FF; % ff correct
end

frames = size(LNCam_filt.i1p848.a425,3);
for k=1:frames
    LNCam_filt.i1p848.a425(:,:,k) = LNCam_filt.i1p848.a425(:,:,k)./LNCam_FF; % ff correct
end

%% Get pixel size

% checkerboard size is 30 mm

%  imtool(mean(LNCam_filt.i1p848.a425, 3),[]);
l = 80;
pixelsize = 30/l;

%% Get cropping rectangle from mean image if needed

LNCam_rect = [600, 500, 399, 399];
RedCam_rect = [600, 500, 399, 399];

%% Crop Images

A_RedCam_filt = struct;
A_LNCam_filt = struct;

frames = size(RedCam_filt.i1p602.a425,3);

A_Red = zeros(400, 400, frames);
A_LN = zeros(400, 400, frames);

for k=1:frames
    A_Red(:,:,k) = imcrop(RedCam_filt.i1p602.a425(:,:,k), RedCam_rect); % crop and assign
    A_LN(:,:,k) = imcrop(LNCam_filt.i1p848.a425(:,:,k), LNCam_rect); % crop and assign
end

A_RedCam_filt.i1p602.a425 = A_Red;  % reassign new array to struct
A_LNCam_filt.i1p848.a425 = A_LN;  % reassign new array to struct


%% NPS Calc Red Cam


m = 400; % size of image
n = 100; % size of each roi
avg_param = 25;

NPS1D_RedCam_filt = struct; % NPS values, will contain a 1D array for each frame of each acquisition
d_RedCam_filt = struct; % mean pixel values, will contain a 1D array for each acquisition stack (1 value per frame)

d_RedCam_filt.i1p602.a425 = zeros(1, frames); % allocate array for frame means (1 value for each frame)
NPS1D_RedCam_filt.i1p602.a425 = zeros(avg_param, frames); % preallocate size for matrix (contains 1D array for each frame)


for k=1:frames

    A_Red   = A_RedCam_filt.i1p602.a425(:,:,k); % currnt frame for nps
    d_RedCam_filt.i1p602.a425(1,k)  = mean(A_Red(:)); % grab mean of frame and assign
    
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
    [NPS_1D, r] = radialavg(NPS, 25);
    NPS_1D = NPS_1D*2; % combine power from positive and negative frequencies

    NPS1D_RedCam_filt.i1p602.a425(:,k) = NPS_1D; % assign array to stuct of means
    
    disp([num2str(100*k/frames), ' % done']);

end


%% NPS Calc LN Cam

m = 400; % size of image
n = 100; % size of each roi
avg_param = 25;

NPS1D_LNCam_filt = struct; % NPS values, will contain a 1D array for each frame of each acquisition
d_LNCam_filt = struct; % mean pixel values, will contain a 1D array for each acquisition stack (1 value per frame)

d_LNCam_filt.i1p848.a425 = zeros(1, frames); % allocate array for frame means (1 value for each frame)
NPS1D_LNCam_filt.i1p848.a425 = zeros(avg_param, frames); % preallocate size for matrix (contains 1D array for each frame)

for k=1:frames
    
    A_LN    = A_LNCam_filt.i1p848.a425(:,:,k);
    d_LNCam_filt.i1p848.a425(1,k)   = mean(A_LN(:)); % grab mean of frame and assign
    
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
    [NPS_1D, r] = radialavg(NPS, 25);
    NPS_1D = NPS_1D*2; % combine power from positive and negative frequencies

    NPS1D_LNCam_filt.i1p848.a425(:,k) = NPS_1D; % assign array to stuct of means
    
    disp([num2str(100*k/frames), ' % done']);

end


%% Average over frames

NPS1D_RedCam_filt_avg = struct;
d_RedCam_filt_avg = struct;

NPS1D_RedCam_filt_avg.i1p602.a425 = mean(NPS1D_RedCam_filt.i1p602.a425,2); % just averaging over frames
d_RedCam_filt_avg.i1p602.a425 = mean(d_RedCam_filt.i1p602.a425,2); % just averaging over frames


NPS1D_LNCam_filt_avg = struct;
d_LNCam_filt_avg = struct;

NPS1D_LNCam_filt_avg.i1p848.a425 = mean(NPS1D_LNCam_filt.i1p848.a425,2); % just averaging over frames
d_LNCam_filt_avg.i1p848.a425 = mean(d_LNCam_filt.i1p848.a425,2); % just averaging over frames


%% Save

save('/Volumes/GoogleDrive/My Drive/Analysis/DQE/DQE Analysis/Results/NPS_filt_results.mat', ...
    'd_RedCam_filt_avg', 'NPS1D_RedCam_filt_avg', 'd_LNCam_filt_avg', 'NPS1D_LNCam_filt_avg', 'r');





