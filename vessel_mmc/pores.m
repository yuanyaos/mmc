% Porous media

clear
% close all

addpath(genpath('/space/neza/2/users/yaoshen/NEU/Research/Redbird/tensorlab/'))
addpath(genpath('/space/neza/2/users/yaoshen/NEU/Research/iso2mesh/'))
addpath('/drives/neza2/users/yaoshen/NEU/Research/mmc/mmcx/mmc/vessel_mmcx/dijkstra')
% [node,face,elem]=meshabox([0 0 0],[50 60 70],100,100);
% 
% prop = 3*ones(size(node,1),1);
% node = [node prop];
% 
% elem = [elem 6*ones(size(elem,1),1) 6*ones(size(elem,1),1) 0*ones(size(elem,1),1) 0*ones(size(elem,1),1)];

%% run mmc

load mesh.mat

clear cfg
cfg.nphoton=1e6;
cfg.node = node;
cfg.elem = elem;
cfg.vessel = 1;
cfg.elemprop=ones(size(cfg.elem,1),1);
cfg.srcpos=[30.5 30.5 0];
cfg.srcdir=[0 0 1];
cfg.prop=[0 0 1 1;0.005 1 0 1.37;0 1 0 1.37];
% cfg.prop=[0 0 1 1;0.036 1.0320 0 1.37;0.21 0.4638 0 1.37];
cfg.tstart=0;
cfg.tend=5e-9;
cfg.tstep=5e-9;
cfg.isreflect=1;
cfg.debuglevel='TP';
% cfg.issaveref=1;  % in addition to volumetric fluence, also save surface diffuse reflectance
cfg.method = 'grid';

flux=mmclab(cfg);
% save flux_vessel_blood flux

figure,slice3(log10(flux.data))
caxis([-5 8])
colormap jet