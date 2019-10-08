% Read radius

clear all
% close all

addpath(genpath('/space/neza/2/users/yaoshen/NEU/Research/Redbird/tensorlab/'))
addpath('/space/neza/2/users/yaoshen/NEU/Research/Optimal_Speration/Tutorial_MMC-MCX/nii')
addpath('/space/neza/2/users/yaoshen/NEU/Research/Optimal_Speration/Tutorial_MMC-MCX/nfri_functions')
addpath('/space/neza/2/users/yaoshen/NEU/Research/Optimal_Speration/Tutorial_MMC-MCX/metch')

vesselr = load_nii('/space/wazu/3/users/yaoshen/backup/iMMC/Vessel-data/radius/99.nii.gz');
r = vesselr.img;

figure,slice3(r(100:200,100:200,100:200))

vesselc = load_nii('/space/wazu/3/users/yaoshen/backup/iMMC/Vessel-data/centerline/99.nii.gz');
center = vesselc.img;

figure,slice3(center(100:200,100:200,100:200))