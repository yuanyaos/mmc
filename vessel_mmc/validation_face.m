% Validation for vessel MMC. Fix the missing face issue

clear
% close all

addpath(genpath('/space/neza/2/users/yaoshen/NEU/Research/Redbird/tensorlab/'))
addpath(genpath('/space/neza/2/users/yaoshen/NEU/Research/iso2mesh/'))
addpath('/drives/neza2/users/yaoshen/NEU/Research/mmc/mmc/vessel_mmc/dijkstra')

%% simple

[node,elem] = meshgrid6(0:1,0:1,0:1);
node = node*60;
node = [node zeros(size(node,1),1)];

figure,
for i=1:size(elem,1)
    plotmesh(node,[elem(i,1) elem(i,2) elem(i,3)],'facealpha',0.1)
    hold on
    plotmesh(node,[elem(i,1) elem(i,2) elem(i,4)],'facealpha',0.1)
    plotmesh(node,[elem(i,1) elem(i,3) elem(i,4)],'facealpha',0.1)
    plotmesh(node,[elem(i,2) elem(i,3) elem(i,4)],'facealpha',0.1)
end
xlabel('x'),ylabel('y'),zlabel('z')

% facet = volface(elem);
surface = [5 8 6; 5 8 7];
[faceelemid, elem_edge] = facelabel(elem,node,surface,0);

faceindex = 6*ones(size(elem,1),1);
faceindex(4) = 1; faceindex(6) = 1;
elem = [elem faceindex 6*ones(size(faceindex,1),3) 3*ones(size(faceindex)) ones(size(faceindex,1),3)];

elem(elem_edge~=0,5) = -elem_edge(elem_edge~=0);
%% run mmc

clear cfg
cfg.nphoton=1e7;
cfg.node = node;
cfg.elem = elem;
cfg.implicit = 2;
cfg.elemprop=ones(size(cfg.elem,1),1);
% cfg.vessel = vessel;
% cfg.radius = 2*ones(size(cfg.elem,1),1);
% cfg.srcpos=[30.5 30.5 60];
cfg.srcpos=[30.5 30.5 60];
cfg.srcdir=[0 0 -1];
cfg.prop=[0 0 1 1;0.005 1 0 1.37;0.005 1 0 6.85];
% cfg.prop=[0 0 1 1;0.005 1 0 1.37;0 0 1 1.37];
% cfg.prop=[0 0 1 1;0.005 1 0 1.37;0.5 1 0 1.37];

cfg.tstart=0;
cfg.tend=5e-9;
cfg.tstep=5e-9;
cfg.isreflect=1;
cfg.debuglevel='TP';
cfg.issaveref=1;  % in addition to volumetric fluence, also save surface diffuse reflectance
cfg.method = 'grid';

cfg.isreflect = 1;

flux=mmclab(cfg);
% save flux_vessel_blood flux

figure,slice3(log10(flux.data))
caxis([-5 8])
colormap jet

figure,imagesc(log10(squeeze(flux.data(:,30,:))),[-5 8]),colormap jet

% save /space/neza/2/users/yaoshen/NEU/Research/mmc/data/face_vessel_n685_r3 flux

%% Compare original MMC and vessel MMC

n = 31;
[xx,yy]=meshgrid(0:59,0:59);
xx = xx'; yy = yy';

levels = -5:8;
figure,
load /drives/neza2/users/yaoshen/NEU/Research/mmc/data/face_original_mu05_r3
flux_original = rot90(squeeze(flux.data(1:end-1,n,1:end-1)),1);
ax1 = subplot(131);
imagesc(log10(flux_original),[-5 8]);
title(['Original, slice=' num2str(n)])
colormap(ax1,jet)

clear flux
load /space/neza/2/users/yaoshen/NEU/Research/mmc/data/face_vessel_mu05_r3
flux_vessel = rot90(squeeze(flux.data(1:end-1,n,1:end-1)),1);
ax2 = subplot(132);
imagesc(log10(flux_vessel),[-5 8]);
title(['Vessel, slice=' num2str(n)])
colormap(ax2,jet)

ax3 = subplot(133);
contourf(xx,yy,rot90(log10(flux_original),-1),levels,'LineWidth',0.1);
hold on
p1 = contour(xx,yy,rot90(log10(flux_original),-1),levels,'-k','ShowText','on');

clear flux
p2 = contour(xx,yy,rot90(log10(flux_vessel),-1),levels,'linestyle','--','color','r','linewidth',1,'ShowText','on');
legend('Original','Original','Vessel')
legend boxoff
% caxis([1.5,5]);
colormap(ax3,parula)