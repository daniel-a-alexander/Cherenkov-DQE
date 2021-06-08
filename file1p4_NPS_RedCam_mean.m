%% Path

addpath(genpath('/Users/danielalexander/OneDrive - Dartmouth College/Matlab Functions'));
datapath = '/Volumes/Extreme1TB/DQE_1p3/';

%% Data Read

% Read in RedCam FF and DF ----------------------------------------------------------------------------------------------------
RedCam_FF                = read_dovi(fullfile(datapath, '2019-09-30 21-52-16-453/meas_s0_cam0.dovi')); 
RedCam_FF                = mean(RedCam_FF, 3);
RedCam_FF_DF             = read_dovi(fullfile(datapath, '2019-09-30 21-52-51-989/meas_s0_cam0.dovi')); % read in DF file for FF
RedCam_FF_DF             = mean(RedCam_FF_DF, 3);
RedCam_FF                = RedCam_FF - RedCam_FF_DF; % need to df first
RedCam_FF(RedCam_FF<0)   = 0;
% RedCam_FF                = medfilt2(RedCam_FF);
RedCam_FF                = RedCam_FF/mean(RedCam_FF(:)); % normalize ff

% Read data ----------------------------------------------------------------------------------------------------

RedCam = struct;

i=0;

RedCam_mean.i1p602.a200 = mean(read_dovi(fullfile(datapath, '2019-09-30 20-41-34-205/meas_s1_cam0.dovi')), 3);i=i+1;disp([num2str(100*i/40), ' % done']);
RedCam_mean.i1p602.a250 = mean(read_dovi(fullfile(datapath, '2019-09-30 20-42-34-658/meas_s1_cam0.dovi')), 3);i=i+1;disp([num2str(100*i/40), ' % done']);
RedCam_mean.i1p602.a300 = mean(read_dovi(fullfile(datapath, '2019-09-30 20-43-41-814/meas_s1_cam0.dovi')), 3);i=i+1;disp([num2str(100*i/40), ' % done']);
RedCam_mean.i1p602.a325 = mean(read_dovi(fullfile(datapath, '2019-09-30 20-44-40-112/meas_s1_cam0.dovi')), 3);i=i+1;disp([num2str(100*i/40), ' % done']);
RedCam_mean.i1p602.a350 = mean(read_dovi(fullfile(datapath, '2019-09-30 20-45-42-385/meas_s1_cam0.dovi')), 3);i=i+1;disp([num2str(100*i/40), ' % done']);
RedCam_mean.i1p602.a375 = mean(read_dovi(fullfile(datapath, '2019-09-30 20-46-57-166/meas_s1_cam0.dovi')), 3);i=i+1;disp([num2str(100*i/40), ' % done']);
RedCam_mean.i1p602.a400 = mean(read_dovi(fullfile(datapath, '2019-09-30 20-47-54-907/meas_s1_cam0.dovi')), 3);i=i+1;disp([num2str(100*i/40), ' % done']);
RedCam_mean.i1p602.a425 = mean(read_dovi(fullfile(datapath, '2019-09-30 20-48-49-269/meas_s1_cam0.dovi')), 3);i=i+1;disp([num2str(100*i/40), ' % done']);

% RedCam_mean.i1p717.a200 = mean(read_dovi(fullfile(datapath, '2019-09-30 21-27-35-173/meas_s1_cam0.dovi')), 3);i=i+1;disp([num2str(100*i/40), ' % done']);
% RedCam_mean.i1p717.a250 = mean(read_dovi(fullfile(datapath, '2019-09-30 21-28-41-396/meas_s1_cam0.dovi')), 3);i=i+1;disp([num2str(100*i/40), ' % done']);
% RedCam_mean.i1p717.a300 = mean(read_dovi(fullfile(datapath, '2019-09-30 21-29-30-620/meas_s1_cam0.dovi')), 3);i=i+1;disp([num2str(100*i/40), ' % done']);
% RedCam_mean.i1p717.a325 = mean(read_dovi(fullfile(datapath, '2019-09-30 21-30-22-798/meas_s1_cam0.dovi')), 3);i=i+1;disp([num2str(100*i/40), ' % done']);
% RedCam_mean.i1p717.a350 = mean(read_dovi(fullfile(datapath, '2019-09-30 21-31-23-049/meas_s1_cam0.dovi')), 3);i=i+1;disp([num2str(100*i/40), ' % done']);
% RedCam_mean.i1p717.a375 = mean(read_dovi(fullfile(datapath, '2019-09-30 21-32-17-845/meas_s1_cam0.dovi')), 3);i=i+1;disp([num2str(100*i/40), ' % done']);
% RedCam_mean.i1p717.a400 = mean(read_dovi(fullfile(datapath, '2019-09-30 21-33-32-754/meas_s1_cam0.dovi')), 3);i=i+1;disp([num2str(100*i/40), ' % done']);
RedCam_mean.i1p717.a425 = mean(read_dovi(fullfile(datapath, '2019-09-30 21-34-32-216/meas_s1_cam0.dovi')), 3);i=i+1;disp([num2str(100*i/40), ' % done']);

% RedCam_mean.i1p847.a200 = mean(read_dovi(fullfile(datapath, '2019-09-30 21-16-32-799/meas_s1_cam0.dovi')), 3);i=i+1;disp([num2str(100*i/40), ' % done']);
% RedCam_mean.i1p847.a250 = mean(read_dovi(fullfile(datapath, '2019-09-30 21-17-20-832/meas_s1_cam0.dovi')), 3);i=i+1;disp([num2str(100*i/40), ' % done']);
% RedCam_mean.i1p847.a300 = mean(read_dovi(fullfile(datapath, '2019-09-30 21-18-13-056/meas_s1_cam0.dovi')), 3);i=i+1;disp([num2str(100*i/40), ' % done']);
% RedCam_mean.i1p847.a325 = mean(read_dovi(fullfile(datapath, '2019-09-30 21-19-13-234/meas_s1_cam0.dovi')), 3);i=i+1;disp([num2str(100*i/40), ' % done']);
% RedCam_mean.i1p847.a350 = mean(read_dovi(fullfile(datapath, '2019-09-30 21-20-40-138/meas_s1_cam0.dovi')), 3);i=i+1;disp([num2str(100*i/40), ' % done']);
% RedCam_mean.i1p847.a375 = mean(read_dovi(fullfile(datapath, '2019-09-30 21-21-35-057/meas_s1_cam0.dovi')), 3);i=i+1;disp([num2str(100*i/40), ' % done']);
% RedCam_mean.i1p847.a400 = mean(read_dovi(fullfile(datapath, '2019-09-30 21-22-36-864/meas_s1_cam0.dovi')), 3);i=i+1;disp([num2str(100*i/40), ' % done']);
RedCam_mean.i1p847.a425 = mean(read_dovi(fullfile(datapath, '2019-09-30 21-23-31-159/meas_s1_cam0.dovi')), 3);i=i+1;disp([num2str(100*i/40), ' % done']);

% RedCam_mean.i1p923.a200 = mean(read_dovi(fullfile(datapath, '2019-09-30 21-03-52-607/meas_s1_cam0.dovi')), 3);i=i+1;disp([num2str(100*i/40), ' % done']);
% RedCam_mean.i1p923.a250 = mean(read_dovi(fullfile(datapath, '2019-09-30 21-05-04-036/meas_s1_cam0.dovi')), 3);i=i+1;disp([num2str(100*i/40), ' % done']);
% RedCam_mean.i1p923.a300 = mean(read_dovi(fullfile(datapath, '2019-09-30 21-06-33-689/meas_s1_cam0.dovi')), 3);i=i+1;disp([num2str(100*i/40), ' % done']);
% RedCam_mean.i1p923.a325 = mean(read_dovi(fullfile(datapath, '2019-09-30 21-08-10-191/meas_s1_cam0.dovi')), 3);i=i+1;disp([num2str(100*i/40), ' % done']);
% RedCam_mean.i1p923.a350 = mean(read_dovi(fullfile(datapath, '2019-09-30 21-09-04-096/meas_s1_cam0.dovi')), 3);i=i+1;disp([num2str(100*i/40), ' % done']);
% RedCam_mean.i1p923.a375 = mean(read_dovi(fullfile(datapath, '2019-09-30 21-10-05-659/meas_s1_cam0.dovi')), 3);i=i+1;disp([num2str(100*i/40), ' % done']);
% RedCam_mean.i1p923.a400 = mean(read_dovi(fullfile(datapath, '2019-09-30 21-11-22-502/meas_s1_cam0.dovi')), 3);i=i+1;disp([num2str(100*i/40), ' % done']);
RedCam_mean.i1p923.a425 = mean(read_dovi(fullfile(datapath, '2019-09-30 21-12-15-194/meas_s1_cam0.dovi')), 3);i=i+1;disp([num2str(100*i/40), ' % done']);

% RedCam_mean.i2p016.a200 = mean(read_dovi(fullfile(datapath, '2019-09-30 20-52-45-925/meas_s1_cam0.dovi')), 3);i=i+1;disp([num2str(100*i/40), ' % done']);
% RedCam_mean.i2p016.a250 = mean(read_dovi(fullfile(datapath, '2019-09-30 20-54-07-093/meas_s1_cam0.dovi')), 3);i=i+1;disp([num2str(100*i/40), ' % done']);
% RedCam_mean.i2p016.a300 = mean(read_dovi(fullfile(datapath, '2019-09-30 20-55-08-226/meas_s1_cam0.dovi')), 3);i=i+1;disp([num2str(100*i/40), ' % done']);
% RedCam_mean.i2p016.a325 = mean(read_dovi(fullfile(datapath, '2019-09-30 20-56-10-659/meas_s1_cam0.dovi')), 3);i=i+1;disp([num2str(100*i/40), ' % done']);
% RedCam_mean.i2p016.a350 = mean(read_dovi(fullfile(datapath, '2019-09-30 20-57-02-236/meas_s1_cam0.dovi')), 3);i=i+1;disp([num2str(100*i/40), ' % done']);
% RedCam_mean.i2p016.a375 = mean(read_dovi(fullfile(datapath, '2019-09-30 20-57-55-971/meas_s1_cam0.dovi')), 3);i=i+1;disp([num2str(100*i/40), ' % done']);
% RedCam_mean.i2p016.a400 = mean(read_dovi(fullfile(datapath, '2019-09-30 20-58-54-916/meas_s1_cam0.dovi')), 3);i=i+1;disp([num2str(100*i/40), ' % done']);
RedCam_mean.i2p016.a425 = mean(read_dovi(fullfile(datapath, '2019-09-30 20-59-55-080/meas_s1_cam0.dovi')), 3);i=i+1;disp([num2str(100*i/40), ' % done']);

clear i

% FF Red Cam averages ----------------------------------------------------------------------------------------------------
intGains_Red = fieldnames(RedCam_mean);

tic
for l=1:numel(intGains_Red) % index l loops over all intensifier gains
    
    anaGains = fieldnames(RedCam_mean.(intGains_Red{l})); % defining this here in case its different for each int gain, shouldnt be
    
    for m=1:numel(anaGains) % index m loops over all analog gains
        
        currentIm = double(RedCam_mean.(intGains_Red{l}).(anaGains{m})); % grab currect acquisition
        FFCIm = currentIm./RedCam_FF;
        RedCam_mean.(intGains_Red{l}).(anaGains{m}) = FFCIm;  % assign FF corrected image to struct
    end
    toc 
end


%% Get pixel size

% checkerboard size is 30 mm

%  imtool(mean(RedCam.i1p848.a425, 3),[]);
l = 80;
pixelsize = 30/l;

%% Get cropping rectangle from mean image if needed

RedCam_rect = [600, 500, 399, 399];

%% Loop through struct to crop

A_RedCam_mean = struct; %cropped images

intGains = fieldnames(RedCam_mean);

tic
for i=1:numel(intGains) % index l loops over all intensifier gains
    
    anaGains = fieldnames(RedCam_mean.(intGains{i})); % defining this here in case its different for each int gain, shouldnt be
    
    for j=1:numel(anaGains) % index m loops over all analog gains
        
        currentIm = double(RedCam_mean.(intGains{i}).(anaGains{j}));
        A_RedCam_mean.(intGains{i}).(anaGains{j}) = imcrop(currentIm, RedCam_rect);  % reassign new array to struct
        
    end
    toc
end

%% NPS calc

m = 400; % size of image
n = 100; % size of each roi
avg_param = 25;

NPS1D_RedCam_mean = struct; % NPS values, will contain a 1D array for each frame of each acquisition
d_RedCam_mean = struct; % mean pixel values, will contain a 1D array for each acquisition stack (1 value per frame)
intGains = fieldnames(A_RedCam_mean);

tic
for i=1:numel(intGains) % index l loops over all intensifier gains
    
    anaGains = fieldnames(A_RedCam_mean.(intGains{i})); % defining this here in case its different for each int gain, shouldnt be
    
    for j=1:numel(anaGains) % index m loops over all analog gains
        
        A = A_RedCam_mean.(intGains{i}).(anaGains{j}); % current image
        
        d_RedCam_mean.(intGains{i}).(anaGains{j}) = mean(A(:)); % allocate array for frame means (1 value for each frame)
        
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

        NPS1D_RedCam_mean.(intGains{i}).(anaGains{j}) = NPS_1D'; % assign array to stuct 

        toc
    end
end

%% Save

save('/Volumes/GoogleDrive/My Drive/Analysis/DQE/DQE Analysis/Results/NPS_RedCam_mean_results.mat', ...
    'd_RedCam_mean', 'NPS1D_RedCam_mean', 'r');


