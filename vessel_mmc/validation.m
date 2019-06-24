clear
close all

addpath(genpath('/space/neza/2/users/yaoshen/NEU/Research/Redbird/tensorlab/'))

[node,elem] = meshgrid6(0:1,0:1,0:1);
node = node*60;

figure,plotmesh(node,elem)

pse(1) = 1;
pse(2) = 8;
[vessel] = vessellabel(elem,node,pse)

%% run mmc

clear cfg
cfg.nphoton=1e7;
cfg.node = node;
cfg.elem = elem;
cfg.elemprop=ones(size(cfg.elem,1),1);
cfg.vessel = 3*ones(size(cfg.elem,1),1);
cfg.radius = 3*ones(size(cfg.elem,1),1);
cfg.srcpos=[35 25 0];
cfg.srcdir=[0 0 1];
cfg.prop=[0 0 1 1;0.005 1 0 1.37;0.005 1 0 13.7];
cfg.tstart=0;
cfg.tend=5e-9;
cfg.tstep=5e-9;
cfg.isreflect=1;
cfg.debuglevel='TP';
cfg.issaveref=1;  % in addition to volumetric fluence, also save surface diffuse reflectance
cfg.method = 'grid';

flux=mmclab(cfg);

figure,slice3(log10(flux.data))
colormap jet