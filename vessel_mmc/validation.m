clear
close all

addpath(genpath('/space/neza/2/users/yaoshen/NEU/Research/Redbird/tensorlab/'))
addpath(genpath('/space/neza/2/users/yaoshen/NEU/Research/iso2mesh/'))
addpath('/drives/neza2/users/yaoshen/NEU/Research/mmc/mmc/vessel_mmc/dijkstra')

% [node,elem] = meshgrid6(0:1,0:1,0:1);
% node = node*60;
% pse(1) = 1;
% pse(2) = 8;
% [vessel] = vessellabel(elem,node,pse);

[node,face,elem]=meshabox([0 0 0],[50 60 70],1000,1000);
elem = elem(:,1:4);
% eleuniq = unique(elem(:));
% pse = round(rand(2,1)*length(eleuniq));
ps(1) = 56;
pe(2) = 2;
pstart = node(pse(1),:);
pend = node(pse(2),:);
[vessel] = vessellabel(elem,node,pse);

figure,plotmesh(node,elem,'facealpha',0.1,'edgealpha',0.5,'facecolor','w')
hold on
plot3(pstart(1,1),pstart(1,2),pstart(1,3),'-og','LineWidth',2)
plot3(pend(end,1),pend(end,2),pend(end,3),'-or','LineWidth',2)
xlabel('x'),ylabel('y'),zlabel('z')
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
% cfg.prop=[0 0 1 1;0.036 1.0320 0 1.37;0.21 0.4638 0 1.37];
cfg.tstart=0;
cfg.tend=5e-9;
cfg.tstep=5e-9;
cfg.isreflect=1;
cfg.debuglevel='TP';
cfg.issaveref=1;  % in addition to volumetric fluence, also save surface diffuse reflectance
cfg.method = 'grid';

flux=mmclab(cfg);
% save flux_vessel_blood flux

figure,slice3(log10(flux.data))
caxis([-5 8])
colormap jet

%% Compare

n = 20;
[xx,yy]=meshgrid(0:60,0:60); 
xx = xx'; yy = yy';

levels = -5:0.5:8;
figure,
load /drives/neza2/users/yaoshen/NEU/Research/mmc/flux_originaln685
flux_original = rot90(squeeze(flux.data(:,:,n)),-1);
ax1 = subplot(131);
imagesc(log10(flux_original),[2 5]);
title(['Original, slice=' num2str(n)])
colormap(ax1,jet)

clear flux
load /drives/neza2/users/yaoshen/NEU/Research/mmc/flux_vesseln685
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