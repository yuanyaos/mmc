%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MMCLAB - Mesh-based Monte Carlo for MATLAB/Octave by Qianqina Fang
%
% In this example, we show the most basic usage of MMCLAB.
%
% This file is part of Mesh-based Monte Carlo (MMC) URL:http://mcx.sf.net/mmc
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% prepare simulation input
addpath('/space/neza/2/users/yaoshen/NEU/Research/iso2mesh/')
addpath(genpath('/space/neza/2/users/yaoshen/NEU/Research/Redbird/tensorlab/'))

clear cfg
cfg.nphoton=1e7;
% [cfg.node, face, cfg.elem]=meshabox([0 0 0],[10 10 10],6);
[cfg.node,cfg.elem] = meshgrid6(0:1,0:1,0:1);
cfg.node = cfg.node*60;
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

%% run the simulation

flux=mmclab(cfg);

%% plotting the result

% plot the cross-section of the fluence
figure,slice3(log10(flux.data))
figure,imagesc((log10(squeeze(flux.data(:,30,:)))))

figure
subplot(121);
plotmesh([cfg.node(:,1:3),log10(abs(flux.data(1:size(cfg.node,1))))],cfg.elem,'y=5','facecolor','interp','linestyle','none')
view([0 1 0]);
colorbar;
caxis([5 7])
subplot(122);
plotmesh(cfg.node,cfg.elem,'y=5')
view([0 1 0]);

% plot the surface diffuse reflectance
if(isfield(cfg,'issaveref') && cfg.issaveref==1)
    subplot(122);
    faces=faceneighbors(cfg.elem,'rowmajor');
    hs=plotmesh(cfg.node,faces,'cdata',log10(flux.dref(:,1)),'linestyle','none');
    colorbar;
end