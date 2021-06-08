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

% Read data ----------------------------------------------------------------------------------------------------

LNCam = struct;

i=0;

LNCam.i1p848.a200 = read_dovi(fullfile(datapath, '2019-09-30 19-46-05-096/meas_s1_cam0.dovi'));i=i+1;disp([num2str(100*i/40), ' % done']);
LNCam.i1p848.a250 = read_dovi(fullfile(datapath, '2019-09-30 19-47-01-704/meas_s1_cam0.dovi'));i=i+1;disp([num2str(100*i/40), ' % done']);
LNCam.i1p848.a300 = read_dovi(fullfile(datapath, '2019-09-30 19-48-58-847/meas_s1_cam0.dovi'));i=i+1;disp([num2str(100*i/40), ' % done']);
LNCam.i1p848.a325 = read_dovi(fullfile(datapath, '2019-09-30 19-49-53-363/meas_s1_cam0.dovi'));i=i+1;disp([num2str(100*i/40), ' % done']);
LNCam.i1p848.a350 = read_dovi(fullfile(datapath, '2019-09-30 19-50-48-072/meas_s1_cam0.dovi'));i=i+1;disp([num2str(100*i/40), ' % done']);
LNCam.i1p848.a375 = read_dovi(fullfile(datapath, '2019-09-30 19-51-46-705/meas_s1_cam0.dovi'));i=i+1;disp([num2str(100*i/40), ' % done']);
LNCam.i1p848.a400 = read_dovi(fullfile(datapath, '2019-09-30 19-52-50-659/meas_s1_cam0.dovi'));i=i+1;disp([num2str(100*i/40), ' % done']);
LNCam.i1p848.a425 = read_dovi(fullfile(datapath, '2019-09-30 19-53-44-158/meas_s1_cam0.dovi'));i=i+1;disp([num2str(100*i/40), ' % done']);

% LNCam.i1p922.a200 = read_dovi(fullfile(datapath, '2019-09-30 20-20-40-625/meas_s1_cam0.dovi'));i=i+1;disp([num2str(100*i/40), ' % done']);
% LNCam.i1p922.a250 = read_dovi(fullfile(datapath, '2019-09-30 20-21-39-997/meas_s1_cam0.dovi'));i=i+1;disp([num2str(100*i/40), ' % done']);
% LNCam.i1p922.a300 = read_dovi(fullfile(datapath, '2019-09-30 20-22-39-729/meas_s1_cam0.dovi'));i=i+1;disp([num2str(100*i/40), ' % done']);
% LNCam.i1p922.a325 = read_dovi(fullfile(datapath, '2019-09-30 20-23-31-105/meas_s1_cam0.dovi'));i=i+1;disp([num2str(100*i/40), ' % done']);
% LNCam.i1p922.a350 = read_dovi(fullfile(datapath, '2019-09-30 20-24-23-483/meas_s1_cam0.dovi'));i=i+1;disp([num2str(100*i/40), ' % done']);
% LNCam.i1p922.a375 = read_dovi(fullfile(datapath, '2019-09-30 20-25-17-361/meas_s1_cam0.dovi'));i=i+1;disp([num2str(100*i/40), ' % done']);
% LNCam.i1p922.a400 = read_dovi(fullfile(datapath, '2019-09-30 20-26-16-297/meas_s1_cam0.dovi'));i=i+1;disp([num2str(100*i/40), ' % done']);
LNCam.i1p922.a425 = read_dovi(fullfile(datapath, '2019-09-30 20-27-18-295/meas_s1_cam0.dovi'));i=i+1;disp([num2str(100*i/40), ' % done']);

% LNCam.i2p034.a200 = read_dovi(fullfile(datapath, '2019-09-30 19-32-57-582/meas_s1_cam0.dovi'));i=i+1;disp([num2str(100*i/40), ' % done']);
% LNCam.i2p034.a250 = read_dovi(fullfile(datapath, '2019-09-30 19-34-04-233/meas_s1_cam0.dovi'));i=i+1;disp([num2str(100*i/40), ' % done']);
% LNCam.i2p034.a300 = read_dovi(fullfile(datapath, '2019-09-30 19-35-27-548/meas_s1_cam0.dovi'));i=i+1;disp([num2str(100*i/40), ' % done']);
% LNCam.i2p034.a325 = read_dovi(fullfile(datapath, '2019-09-30 19-36-20-311/meas_s1_cam0.dovi'));i=i+1;disp([num2str(100*i/40), ' % done']);
% LNCam.i2p034.a350 = read_dovi(fullfile(datapath, '2019-09-30 19-37-17-848/meas_s1_cam0.dovi'));i=i+1;disp([num2str(100*i/40), ' % done']);
% LNCam.i2p034.a375 = read_dovi(fullfile(datapath, '2019-09-30 19-38-11-963/meas_s1_cam0.dovi'));i=i+1;disp([num2str(100*i/40), ' % done']);
% LNCam.i2p034.a400 = read_dovi(fullfile(datapath, '2019-09-30 19-39-50-995/meas_s1_cam0.dovi'));i=i+1;disp([num2str(100*i/40), ' % done']);
LNCam.i2p034.a425 = read_dovi(fullfile(datapath, '2019-09-30 19-41-22-846/meas_s1_cam0.dovi'));i=i+1;disp([num2str(100*i/40), ' % done']);

% LNCam.i2p135.a200 = read_dovi(fullfile(datapath, '2019-09-30 19-59-54-593/meas_s1_cam0.dovi'));i=i+1;disp([num2str(100*i/40), ' % done']);
% LNCam.i2p135.a250 = read_dovi(fullfile(datapath, '2019-09-30 20-00-52-637/meas_s1_cam0.dovi'));i=i+1;disp([num2str(100*i/40), ' % done']);
% LNCam.i2p135.a300 = read_dovi(fullfile(datapath, '2019-09-30 20-01-54-164/meas_s1_cam0.dovi'));i=i+1;disp([num2str(100*i/40), ' % done']);
% LNCam.i2p135.a325 = read_dovi(fullfile(datapath, '2019-09-30 20-02-45-878/meas_s1_cam0.dovi'));i=i+1;disp([num2str(100*i/40), ' % done']);
% LNCam.i2p135.a350 = read_dovi(fullfile(datapath, '2019-09-30 20-03-45-502/meas_s1_cam0.dovi'));i=i+1;disp([num2str(100*i/40), ' % done']);
% LNCam.i2p135.a375 = read_dovi(fullfile(datapath, '2019-09-30 20-04-54-255/meas_s1_cam0.dovi'));i=i+1;disp([num2str(100*i/40), ' % done']);
% LNCam.i2p135.a400 = read_dovi(fullfile(datapath, '2019-09-30 20-06-23-200/meas_s1_cam0.dovi'));i=i+1;disp([num2str(100*i/40), ' % done']);
LNCam.i2p135.a425 = read_dovi(fullfile(datapath, '2019-09-30 20-07-13-294/meas_s1_cam0.dovi'));i=i+1;disp([num2str(100*i/40), ' % done']);

% LNCam.i2p256.a200 = read_dovi(fullfile(datapath, '2019-09-30 20-09-27-580/meas_s1_cam0.dovi'));i=i+1;disp([num2str(100*i/40), ' % done']);
% LNCam.i2p256.a250 = read_dovi(fullfile(datapath, '2019-09-30 20-10-25-908/meas_s1_cam0.dovi'));i=i+1;disp([num2str(100*i/40), ' % done']);
% LNCam.i2p256.a300 = read_dovi(fullfile(datapath, '2019-09-30 20-11-34-236/meas_s1_cam0.dovi'));i=i+1;disp([num2str(100*i/40), ' % done']);
% LNCam.i2p256.a325 = read_dovi(fullfile(datapath, '2019-09-30 20-12-53-438/meas_s1_cam0.dovi'));i=i+1;disp([num2str(100*i/40), ' % done']);
% LNCam.i2p256.a350 = read_dovi(fullfile(datapath, '2019-09-30 20-14-00-377/meas_s1_cam0.dovi'));i=i+1;disp([num2str(100*i/40), ' % done']);
% LNCam.i2p256.a375 = read_dovi(fullfile(datapath, '2019-09-30 20-15-08-868/meas_s1_cam0.dovi'));i=i+1;disp([num2str(100*i/40), ' % done']);
% LNCam.i2p256.a400 = read_dovi(fullfile(datapath, '2019-09-30 20-15-59-840/meas_s1_cam0.dovi'));i=i+1;disp([num2str(100*i/40), ' % done']);
LNCam.i2p256.a425 = read_dovi(fullfile(datapath, '2019-09-30 20-16-50-844/meas_s1_cam0.dovi'));i=i+1;disp([num2str(100*i/40), ' % done']);

clear i

%% FF LN Cam raw stacks ----------------------------------------------------------------------------------------------------
intGains_LN = fieldnames(LNCam);

tic
for l=1:numel(intGains_LN) % index l loops over all intensifier gains
    
    anaGains = fieldnames(LNCam.(intGains_LN{l})); % defining this here in case its different for each int gain, shouldnt be
    
    for m=1:numel(anaGains) % index m loops over all analog gains
        currentStack = double(LNCam.(intGains_LN{l}).(anaGains{m})); % grab currect acquisition
        ffcStack = zeros(size(currentStack));
        frames = size(currentStack,3);

        for k=1:frames
            ffcStack(:,:,k)  = currentStack(:,:,k)./LNCam_FF; % ff correct
%             temp(~isfinite(temp)) = nan;
%             temp = inpaint_nans(temp);
%             ffcStack(:,:,k) = temp;
            if mod(round(100*k/frames), 10) ==0
                disp(['int gain ' , intGains_LN{l},', cmos gain ', anaGains{m}, ': ' ,num2str(round(100*k/frames)), ' % done']);
            end
        end
        LNCam.(intGains_LN{l}).(anaGains{m}) = ffcStack;  % reassign new array to struct

    end

end

clear m k l anaGains intGains_LN frames

%% Get pixel size

% checkerboard size is 30 mm

%  imtool(mean(LNCam.i1p848.a425, 3),[]);
l = 80;
pixelsize = 30/l;

%% Get cropping rectangle from mean image if needed

LNCam_rect = [600, 500, 399, 399];

%% Loop through struct to crop

A_LNCam = struct; %cropped images

intGains = fieldnames(LNCam);

tic
for i=1:numel(intGains) % index l loops over all intensifier gains
    
    anaGains = fieldnames(LNCam.(intGains{i})); % defining this here in case its different for each int gain, shouldnt be
    
    for j=1:numel(anaGains) % index m loops over all analog gains
        
        currentStack = double(LNCam.(intGains{i}).(anaGains{j})); % grab currect acquisition
        frames = size(currentStack,3);
        A = zeros(400, 400, frames);
        
        for k=1:frames
            A(:,:,k) = imcrop(currentStack(:,:,k), LNCam_rect); % crop and assign
        end
        
        A_LNCam.(intGains{i}).(anaGains{j}) = A;  % reassign new array to struct
        
    end
    toc
end

%% Clear original stack

clear LNCam

%% NPS calc

m = 400; % size of image
n = 100; % size of each roi
avg_param = 25;

NPS1D_LNCam = struct; % NPS values, will contain a 1D array for each frame of each acquisition
d_LNCam = struct; % mean pixel values, will contain a 1D array for each acquisition stack (1 value per frame)
intGains = fieldnames(A_LNCam);

for i=1:numel(intGains) % index l loops over all intensifier gains
    
    anaGains = fieldnames(A_LNCam.(intGains{i})); % defining this here in case its different for each int gain, shouldnt be
    
    for j=1:numel(anaGains) % index m loops over all analog gains
        
        A_stack = A_LNCam.(intGains{i}).(anaGains{j}); % current stack of cropped images, length is number of frames
        frames = size(A_stack,3);
        
        d_LNCam.(intGains{i}).(anaGains{j}) = zeros(1, frames); % allocate array for frame means (1 value for each frame)
        NPS1D_LNCam.(intGains{i}).(anaGains{j}) = zeros(avg_param, frames); % preallocate size for matrix (contains 1D array for each frame)
        
        for k=1:frames
            
            A = A_stack(:,:,k); % currnt frame for nps
            
            d_LNCam.(intGains{i}).(anaGains{j})(1,k) = mean(A(:)); % grab mean of frame and assign
            
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
            
            NPS1D_LNCam.(intGains{i}).(anaGains{j})(:,k) = NPS_1D; % assign array to stuct of means (this is all we're saving here, other than the mean d)
            
            disp(['int gain ' , intGains{i},', cmos gain ', anaGains{j}, ': ' ,num2str(100*k/frames), ' % done']);
        end
    end
end

%% Average out all 1D NPS measurements over frames

NPS1D_LNCam_avg = struct; % Average 1D NPS, will contain a 1D array for each acquisition
d_LNCam_avg = struct; % mean pixel values, will contain a value for each acquisition 

for i=1:numel(intGains) % index l loops over all intensifier gains
    
    anaGains = fieldnames(A_LNCam.(intGains{i})); % defining this here in case its different for each int gain, shouldnt be
    
    for j=1:numel(anaGains) % index m loops over all analog gains
        
        NPS1D_LNCam_avg.(intGains{i}).(anaGains{j}) = mean(NPS1D_LNCam.(intGains{i}).(anaGains{j}),2); % just averaging over frames
        d_LNCam_avg.(intGains{i}).(anaGains{j}) = mean(d_LNCam.(intGains{i}).(anaGains{j}),2); % just averaging over frames
        
    end
end

%% Save

save('/Volumes/GoogleDrive/My Drive/Analysis/DQE/DQE Analysis/Results/NPS_LNCam_results.mat', ...
    'd_LNCam_avg', 'NPS1D_LNCam_avg', 'r');


