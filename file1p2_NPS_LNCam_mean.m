%% Path

addpath(genpath('/Volumes/Google Drove/My Drive/Matlab Functions/'));
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

% Read data ----------------------------------------------------------------------------------------------------

LNCam_mean = struct;

i=0;

LNCam_mean.i1p848.a200 = mean(read_dovi(fullfile(datapath, '2019-09-30 19-46-05-096/meas_s1_cam0.dovi')), 3);i=i+1;disp([num2str(100*i/40), ' % done']);
LNCam_mean.i1p848.a250 = mean(read_dovi(fullfile(datapath, '2019-09-30 19-47-01-704/meas_s1_cam0.dovi')), 3);i=i+1;disp([num2str(100*i/40), ' % done']);
LNCam_mean.i1p848.a300 = mean(read_dovi(fullfile(datapath, '2019-09-30 19-48-58-847/meas_s1_cam0.dovi')), 3);i=i+1;disp([num2str(100*i/40), ' % done']);
LNCam_mean.i1p848.a325 = mean(read_dovi(fullfile(datapath, '2019-09-30 19-49-53-363/meas_s1_cam0.dovi')), 3);i=i+1;disp([num2str(100*i/40), ' % done']);
LNCam_mean.i1p848.a350 = mean(read_dovi(fullfile(datapath, '2019-09-30 19-50-48-072/meas_s1_cam0.dovi')), 3);i=i+1;disp([num2str(100*i/40), ' % done']);
LNCam_mean.i1p848.a375 = mean(read_dovi(fullfile(datapath, '2019-09-30 19-51-46-705/meas_s1_cam0.dovi')), 3);i=i+1;disp([num2str(100*i/40), ' % done']);
LNCam_mean.i1p848.a400 = mean(read_dovi(fullfile(datapath, '2019-09-30 19-52-50-659/meas_s1_cam0.dovi')), 3);i=i+1;disp([num2str(100*i/40), ' % done']);
LNCam_mean.i1p848.a425 = mean(read_dovi(fullfile(datapath, '2019-09-30 19-53-44-158/meas_s1_cam0.dovi')), 3);i=i+1;disp([num2str(100*i/40), ' % done']);

LNCam_mean.i1p922.a200 = mean(read_dovi(fullfile(datapath, '2019-09-30 20-20-40-625/meas_s1_cam0.dovi')), 3);i=i+1;disp([num2str(100*i/40), ' % done']);
LNCam_mean.i1p922.a250 = mean(read_dovi(fullfile(datapath, '2019-09-30 20-21-39-997/meas_s1_cam0.dovi')), 3);i=i+1;disp([num2str(100*i/40), ' % done']);
LNCam_mean.i1p922.a300 = mean(read_dovi(fullfile(datapath, '2019-09-30 20-22-39-729/meas_s1_cam0.dovi')), 3);i=i+1;disp([num2str(100*i/40), ' % done']);
LNCam_mean.i1p922.a325 = mean(read_dovi(fullfile(datapath, '2019-09-30 20-23-31-105/meas_s1_cam0.dovi')), 3);i=i+1;disp([num2str(100*i/40), ' % done']);
LNCam_mean.i1p922.a350 = mean(read_dovi(fullfile(datapath, '2019-09-30 20-24-23-483/meas_s1_cam0.dovi')), 3);i=i+1;disp([num2str(100*i/40), ' % done']);
LNCam_mean.i1p922.a375 = mean(read_dovi(fullfile(datapath, '2019-09-30 20-25-17-361/meas_s1_cam0.dovi')), 3);i=i+1;disp([num2str(100*i/40), ' % done']);
LNCam_mean.i1p922.a400 = mean(read_dovi(fullfile(datapath, '2019-09-30 20-26-16-297/meas_s1_cam0.dovi')), 3);i=i+1;disp([num2str(100*i/40), ' % done']);
LNCam_mean.i1p922.a425 = mean(read_dovi(fullfile(datapath, '2019-09-30 20-27-18-295/meas_s1_cam0.dovi')), 3);i=i+1;disp([num2str(100*i/40), ' % done']);

LNCam_mean.i2p034.a200 = mean(read_dovi(fullfile(datapath, '2019-09-30 19-32-57-582/meas_s1_cam0.dovi')), 3);i=i+1;disp([num2str(100*i/40), ' % done']);
LNCam_mean.i2p034.a250 = mean(read_dovi(fullfile(datapath, '2019-09-30 19-34-04-233/meas_s1_cam0.dovi')), 3);i=i+1;disp([num2str(100*i/40), ' % done']);
LNCam_mean.i2p034.a300 = mean(read_dovi(fullfile(datapath, '2019-09-30 19-35-27-548/meas_s1_cam0.dovi')), 3);i=i+1;disp([num2str(100*i/40), ' % done']);
LNCam_mean.i2p034.a325 = mean(read_dovi(fullfile(datapath, '2019-09-30 19-36-20-311/meas_s1_cam0.dovi')), 3);i=i+1;disp([num2str(100*i/40), ' % done']);
LNCam_mean.i2p034.a350 = mean(read_dovi(fullfile(datapath, '2019-09-30 19-37-17-848/meas_s1_cam0.dovi')), 3);i=i+1;disp([num2str(100*i/40), ' % done']);
LNCam_mean.i2p034.a375 = mean(read_dovi(fullfile(datapath, '2019-09-30 19-38-11-963/meas_s1_cam0.dovi')), 3);i=i+1;disp([num2str(100*i/40), ' % done']);
LNCam_mean.i2p034.a400 = mean(read_dovi(fullfile(datapath, '2019-09-30 19-39-50-995/meas_s1_cam0.dovi')), 3);i=i+1;disp([num2str(100*i/40), ' % done']);
LNCam_mean.i2p034.a425 = mean(read_dovi(fullfile(datapath, '2019-09-30 19-41-22-846/meas_s1_cam0.dovi')), 3);i=i+1;disp([num2str(100*i/40), ' % done']);

LNCam_mean.i2p135.a200 = mean(read_dovi(fullfile(datapath, '2019-09-30 19-59-54-593/meas_s1_cam0.dovi')), 3);i=i+1;disp([num2str(100*i/40), ' % done']);
LNCam_mean.i2p135.a250 = mean(read_dovi(fullfile(datapath, '2019-09-30 20-00-52-637/meas_s1_cam0.dovi')), 3);i=i+1;disp([num2str(100*i/40), ' % done']);
LNCam_mean.i2p135.a300 = mean(read_dovi(fullfile(datapath, '2019-09-30 20-01-54-164/meas_s1_cam0.dovi')), 3);i=i+1;disp([num2str(100*i/40), ' % done']);
LNCam_mean.i2p135.a325 = mean(read_dovi(fullfile(datapath, '2019-09-30 20-02-45-878/meas_s1_cam0.dovi')), 3);i=i+1;disp([num2str(100*i/40), ' % done']);
LNCam_mean.i2p135.a350 = mean(read_dovi(fullfile(datapath, '2019-09-30 20-03-45-502/meas_s1_cam0.dovi')), 3);i=i+1;disp([num2str(100*i/40), ' % done']);
LNCam_mean.i2p135.a375 = mean(read_dovi(fullfile(datapath, '2019-09-30 20-04-54-255/meas_s1_cam0.dovi')), 3);i=i+1;disp([num2str(100*i/40), ' % done']);
LNCam_mean.i2p135.a400 = mean(read_dovi(fullfile(datapath, '2019-09-30 20-06-23-200/meas_s1_cam0.dovi')), 3);i=i+1;disp([num2str(100*i/40), ' % done']);
LNCam_mean.i2p135.a425 = mean(read_dovi(fullfile(datapath, '2019-09-30 20-07-13-294/meas_s1_cam0.dovi')), 3);i=i+1;disp([num2str(100*i/40), ' % done']);

LNCam_mean.i2p256.a200 = mean(read_dovi(fullfile(datapath, '2019-09-30 20-09-27-580/meas_s1_cam0.dovi')), 3);i=i+1;disp([num2str(100*i/40), ' % done']);
LNCam_mean.i2p256.a250 = mean(read_dovi(fullfile(datapath, '2019-09-30 20-10-25-908/meas_s1_cam0.dovi')), 3);i=i+1;disp([num2str(100*i/40), ' % done']);
LNCam_mean.i2p256.a300 = mean(read_dovi(fullfile(datapath, '2019-09-30 20-11-34-236/meas_s1_cam0.dovi')), 3);i=i+1;disp([num2str(100*i/40), ' % done']);
LNCam_mean.i2p256.a325 = mean(read_dovi(fullfile(datapath, '2019-09-30 20-12-53-438/meas_s1_cam0.dovi')), 3);i=i+1;disp([num2str(100*i/40), ' % done']);
LNCam_mean.i2p256.a350 = mean(read_dovi(fullfile(datapath, '2019-09-30 20-14-00-377/meas_s1_cam0.dovi')), 3);i=i+1;disp([num2str(100*i/40), ' % done']);
LNCam_mean.i2p256.a375 = mean(read_dovi(fullfile(datapath, '2019-09-30 20-15-08-868/meas_s1_cam0.dovi')), 3);i=i+1;disp([num2str(100*i/40), ' % done']);
LNCam_mean.i2p256.a400 = mean(read_dovi(fullfile(datapath, '2019-09-30 20-15-59-840/meas_s1_cam0.dovi')), 3);i=i+1;disp([num2str(100*i/40), ' % done']);
LNCam_mean.i2p256.a425 = mean(read_dovi(fullfile(datapath, '2019-09-30 20-16-50-844/meas_s1_cam0.dovi')), 3);i=i+1;disp([num2str(100*i/40), ' % done']);

clear i

% FF LN Cam averages ----------------------------------------------------------------------------------------------------
intGains_LN = fieldnames(LNCam_mean);

tic
for l=1:numel(intGains_LN) % index l loops over all intensifier gains
    
    anaGains = fieldnames(LNCam_mean.(intGains_LN{l})); % defining this here in case its different for each int gain, shouldnt be
    
    for m=1:numel(anaGains) % index m loops over all analog gains
        
        currentIm = double(LNCam_mean.(intGains_LN{l}).(anaGains{m})); % grab currect acquisition
        FFCIm = currentIm./LNCam_FF;
        LNCam_mean.(intGains_LN{l}).(anaGains{m}) = FFCIm;  % assign FF corrected image to struct
        
    end

end

%% Get pixel size

% checkerboard size is 30 mm

%  imtool(mean(LNCam.i1p848.a425, 3),[]);
l = 80;
pixelsize = 30/l;

%% Get cropping rectangle from mean image if needed

LNCam_rect = [600, 500, 399, 399];

%% Loop through struct to crop

A_LNCam_mean = struct; %cropped images

intGains = fieldnames(LNCam_mean);

tic
for i=1:numel(intGains) % index l loops over all intensifier gains
    
    anaGains = fieldnames(LNCam_mean.(intGains{i})); % defining this here in case its different for each int gain, shouldnt be
    
    for j=1:numel(anaGains) % index m loops over all analog gains
        
        currentIm = double(LNCam_mean.(intGains{i}).(anaGains{j}));
        A_LNCam_mean.(intGains{i}).(anaGains{j}) = imcrop(currentIm, LNCam_rect);  % reassign new array to struct
        
    end
    toc
end

%% NPS calc

m = 400; % size of image
n = 100; % size of each roi
avg_param = 25;

NPS1D_LNCam_mean = struct; % NPS values, will contain a 1D array for each frame of each acquisition
d_LNCam_mean = struct; % mean pixel values, will contain a 1D array for each acquisition stack (1 value per frame)
intGains = fieldnames(A_LNCam_mean);

tic
for i=1:numel(intGains) % index l loops over all intensifier gains
    
    anaGains = fieldnames(A_LNCam_mean.(intGains{i})); % defining this here in case its different for each int gain, shouldnt be
    
    for j=1:numel(anaGains) % index m loops over all analog gains
        
        A = A_LNCam_mean.(intGains{i}).(anaGains{j}); % current image
        
        d_LNCam_mean.(intGains{i}).(anaGains{j}) = mean(A(:)); % allocate array for frame means (1 value for each frame)
        
        FFTstack = zeros(n,n,length(1:n/2:m-n/2)^2); % this is the FFT for each ROI
    
        z = 1; % ROI index

        for g = 1:n/2:m-n/2
            for h = 1:n/2:m-n/2
                rect_gh = [g, h, n-1, n-1]; 
                Asub = imcrop(A, rect_gh); % get ROI from cropped image 
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

        NPS1D_LNCam_mean.(intGains{i}).(anaGains{j}) = NPS_1D'; % assign array to stuct 

        toc
    end
end

%% Save

save('/Volumes/GoogleDrive/My Drive/Analysis/DQE/DQE Analysis/Results/NPS_LNCam_mean_results.mat', ...
    'd_LNCam_mean', 'NPS1D_LNCam_mean', 'r');



