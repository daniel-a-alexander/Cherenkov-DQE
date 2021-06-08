
%% Add Path

addpath(genpath('/Volumes/GoogleDrive/My Drive/Analysis/DQE/DQE Analysis'));
addpath(genpath('/Volumes/GoogleDrive/My Drive/Matlab Functions'));

load('/Volumes/GoogleDrive/My Drive/Analysis/DQE/DQE Analysis/Results/SpectFluenceData.mat');


%% q from laser

pixelsizelaser  = 0.3552; % from DQE_1p1
Area            = (512*pixelsizelaser)^2; % Area of illumination
Power           = 0.005e-6; % .005 uW
DutyRatio       = 0.14/6; % Ratio of duty cycles
PowerC          = Power*DutyRatio; % corrected power
T               = 51e-3; % 51 ms exposure time
E               = 3.094e-19; % 642 nm photon energy in Joules

q_laser               = PowerC*T/(E*Area);

%% Change QE to fraction

QE(:,2) = QE(:,2)./100;

%% Downsample cherenkov

CherenkovSpectrum2 = zeros(size(QE));

CherenkovSpectrum2(:,1) = QE(:,1);

CherenkovSpectrum2(:,2) = interp1(CherenkovSpectrum(:,1), CherenkovSpectrum(:,2), CherenkovSpectrum2(:,1));

%% Multiply Curves

FluenceSpectrum = zeros(size(QE));
FluenceSpectrum(:,1) = QE(:,1);

FluenceSpectrum(:,2) = QE(:,2).*CherenkovSpectrum2(:,2);

%% Load pixel data

load('/Volumes/GoogleDrive/My Drive/Analysis/DQE/DQE Analysis/Results/d_laser.mat', 'd_1p603V');
load('/Volumes/GoogleDrive/My Drive/Analysis/DQE/DQE Analysis/Results/NPS_all_results.mat', 'd_RedCam_avg');

% Pixel counts at gain 400
LaserCounts = d_1p603V(7);
CherenkovCounts = d_RedCam_avg.i1p602.a400;

clear d_1p603V d_RedCam_avg CherenkovSpectrum

%% Find 642 nm in spectrum

% Normalize spectrum by value at wl closest to 642 (640)
FluenceSpectrum(:,2) = FluenceSpectrum(:,2)./FluenceSpectrum(49,2);

%% Now scale by the measured q/count

FluenceCountRatio  = q_laser/LaserCounts;

FluenceSpectrum(:,2) = FluenceSpectrum(:,2).*FluenceCountRatio;

%% Plot

figure;
plot(FluenceSpectrum(:,1),FluenceSpectrum(:,2));

%%

AverageFluence = trapz(FluenceSpectrum(:,1),FluenceSpectrum(:,2));











