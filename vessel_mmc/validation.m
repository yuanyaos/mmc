clear
close all

addpath(genpath('/space/neza/2/users/yaoshen/NEU/Research/Redbird/tensorlab/'))
addpath(genpath('/space/neza/2/users/yaoshen/NEU/Research/iso2mesh/'))

[node,elem] = meshgrid6(0:1,0:1,0:1);
node = node*60;

% figure,plotmesh(node,elem)

pse(1) = 1;
pse(2) = 8;
[vessel] = vessellabel(elem,node,pse);

%% run mmc

clear cfg
cfg.nphoton=1e7;
cfg.node = node;
cfg.elem = elem;
cfg.elemprop=ones(size(cfg.elem,1),1);
cfg.vessel = vessel;
cfg.radius = 2*ones(size(cfg.elem,1),1);
cfg.srcpos=[30.5 30.5 0];
cfg.srcdir=[0 0 1];
cfg.prop=[0 0 1 1;0.005 1 0 1.37;0.5 1 0 1.37];
cfg.tstart=0;
cfg.tend=5e-9;
cfg.tstep=5e-9;
cfg.isreflect=1;
cfg.debuglevel='TP';
cfg.issaveref=1;  % in addition to volumetric fluence, also save surface diffuse reflectance
cfg.method = 'grid';

flux=mmclab(cfg);
save flux_vessel flux

figure,slice3(log10(flux.data))
caxis([-5 8])
colormap jet

%% Compare

n = 40;
[xx,yy]=meshgrid(0:60,0:60); 
xx = xx'; yy = yy';

levels = -5:0.5:8;
figure,
load /drives/neza2/users/yaoshen/NEU/Research/mmc/flux_original
flux_original = rot90(squeeze(flux.data(:,:,n)),-1);
ax1 = subplot(131);
imagesc(log10(flux_original),[2 5]);
title(['Original, slice=' num2str(n)])
colormap(ax1,jet)

clear flux
load /drives/neza2/users/yaoshen/NEU/Research/mmc/flux_vessel
flux_vessel = rot90(squeeze(flux.data(:,:,n)),-1);
ax2 = subplot(132);
imagesc(log10(flux_vessel),[2 5]);
title(['Vessel, slice=' num2str(n)])
colormap(ax2,jet)

ax3 = subplot(133);
contourf(xx,yy,log10(flux_original),levels,'LineWidth',0.1);
hold on
p1 = contour(xx,yy,log10(flux_original),levels,'-k','ShowText','on');

clear flux
p2 = contour(xx,yy,log10(flux_vessel),levels,'linestyle','--','color','w','linewidth',1,'ShowText','on');
legend('Original','Original','Vessel')
legend boxoff
caxis([1.5,5]);
colormap(ax3,parula)