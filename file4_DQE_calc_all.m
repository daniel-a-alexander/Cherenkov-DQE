%% LOAD DATASET HERE!!!! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

load('/Users/danielalexander/Documents/Analysis/DQE/DQE Analysis/Results/NPS_all_results.mat'); % Loads up all NPS and pixel count data
load('/Users/danielalexander/Documents/Analysis/DQE/DQE Analysis/Results/MTF_results.mat'); % loads up frequency domain and MTF data
load('/Users/danielalexander/Documents/Analysis/DQE/DQE Analysis/Results/q_cherenkov.mat'); % loads up input quanta data

%% Calc DQE Red Cam

DQE_RedCam = struct;

intGains = fieldnames(d_RedCam_avg);

for i=1:numel(intGains) % index l loops over all intensifier gains
    
    anaGains = fieldnames(d_RedCam_avg.(intGains{i}));
    
    for j=1:numel(anaGains) % index m loops over all analog gains
        
        DQE_RedCam.(intGains{i}).(anaGains{j}) = d_RedCam_avg.(intGains{i}).(anaGains{j})^2 * (MTF.RedCam.nofilt.mtf.^2) ...
            ./ ( q_cherenkov * NPS1D_RedCam_avg.(intGains{i}).(anaGains{j}) ); 
        
    end
end

%% Calc DQE Red Cam filt

DQE_RedCam_filt = struct;

%Fluence is multiplied by 5 for 5-frame filter

DQE_RedCam_filt.i1p602.a425 = d_RedCam_filt_avg.i1p602.a425^2 * (MTF.RedCam.filt.mtf.^2) ...
    ./ ( 5*q_cherenkov * NPS1D_RedCam_filt_avg.i1p602.a425 ); % Fluence is multiplied by 5 for 5-frame filter

%% Calc DQE Red Cam Mean

DQE_RedCam_mean = struct;

intGains = fieldnames(d_RedCam_mean);

for i=1:numel(intGains) % index l loops over all intensifier gains
    
    anaGains = fieldnames(d_RedCam_mean.(intGains{i}));
    
    for j=1:numel(anaGains) % index m loops over all analog gains
        
        DQE_RedCam_mean.(intGains{i}).(anaGains{j}) = d_RedCam_mean.(intGains{i}).(anaGains{j})^2 * (MTF.RedCam.nofilt.mtf.^2) ...
            ./ ( 100*q_cherenkov * NPS1D_RedCam_mean.(intGains{i}).(anaGains{j}) ); % Fluence is multiplied by 100 for 100-frame average
        
    end
end




%% Calc DQE LN Cam

DQE_LNCam = struct;

intGains = fieldnames(d_LNCam_avg);

for i=1:numel(intGains) % index l loops over all intensifier gains
    
    anaGains = fieldnames(d_LNCam_avg.(intGains{i}));
    
    for j=1:numel(anaGains) % index m loops over all analog gains
        
        DQE_LNCam.(intGains{i}).(anaGains{j}) = d_LNCam_avg.(intGains{i}).(anaGains{j})^2 * (MTF.LNCam.nofilt.mtf.^2) ...
            ./ ( q_cherenkov * NPS1D_LNCam_avg.(intGains{i}).(anaGains{j}) ); 
        
    end
end


%% Calc DQE LN Cam filt

DQE_LNCam_filt = struct;


DQE_LNCam_filt.i1p848.a425 = d_LNCam_filt_avg.i1p848.a425^2 * (MTF.LNCam.filt.mtf.^2) ...
    ./ ( 5*q_cherenkov * NPS1D_LNCam_filt_avg.i1p848.a425 ); % Fluence is multiplied by 5 for 5-frame filter


%% Calc DQE LN Cam Mean 

DQE_LNCam_mean = struct;

intGains = fieldnames(d_LNCam_mean);

for i=1:numel(intGains) % index l loops over all intensifier gains
    
    anaGains = fieldnames(d_LNCam_mean.(intGains{i}));
    
    for j=1:numel(anaGains) % index m loops over all analog gains
        
        DQE_LNCam_mean.(intGains{i}).(anaGains{j}) = d_LNCam_mean.(intGains{i}).(anaGains{j})^2 * (MTF.LNCam.nofilt.mtf.^2) ...
            ./ ( 100*q_cherenkov * NPS1D_LNCam_mean.(intGains{i}).(anaGains{j}) ); % Fluence is multiplied by 100 for 100-frame average
        
    end
end

%% Clear

clearvars anaGains intGains i j

%% Save all variables

save('/Users/danielalexander/Documents/Analysis/DQE/DQE Analysis/Results/DQE_all_results.mat');
