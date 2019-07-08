% Validation for vessel MMC

clear
% close all

addpath(genpath('/space/neza/2/users/yaoshen/NEU/Research/Redbird/tensorlab/'))
addpath(genpath('/space/neza/2/users/yaoshen/NEU/Research/iso2mesh/'))
addpath('/drives/neza2/users/yaoshen/NEU/Research/mmc/mmc/vessel_mmc/dijkstra')

%% simple

[node,elem] = meshgrid6(0:1,0:1,0:1);
node = node*60;
node = [node [1 0 0 0 0 0 0 1]'];
pse(1) = 1;
pse(2) = 8;
[vessel,vesseln,pathelem,nodeelem] = vessellabel(elem,node,pse,2,0);
elem = [elem vessel 6*ones(size(vessel)) ones(size(vessel)) ones(size(vessel))];

%% complex

% [node,face,elem]=meshabox([0 0 0],[50 60 70],1000,1000);
% elem = elem(:,1:4);
% % eleuniq = unique(elem(:));
% % pse = round(rand(2,1)*length(eleuniq));
% pse(1) = 56;
% pse(2) = 2;
% [vessel1,vesseln1,pathelem1,nodeelem1] = vessellabel(elem,node,pse,2,0);
% 
% pse(1) = 70;
% pse(2) = 8;
% [vessel2,vesseln2,pathelem2,nodeelem2] = vessellabel(elem,node,pse,1,0);
% 
% % combine
% vd = vessel1-vessel2;
% i1 = find(vd==0);
% i2 = find(vessel2~=6);
% i2 = intersect(i1,i2);  % overlap local index
% vessel2(i2) = 6;
% vessel = [vessel1 vessel2 2*ones(size(vessel1)) 1*ones(size(vessel2))];
% 
% vesseln = [vesseln1, vesseln2];
% vesseln = max(vesseln,[],2);

%% run mmc

clear cfg
cfg.nphoton=1e6;
cfg.node = node;
% cfg.node = [cfg.node vesseln];
cfg.elem = elem;
% cfg.elem = [cfg.elem vessel];
% plotvessel(cfg.elem,cfg.node)
cfg.vessel = 1;
cfg.elemprop=ones(size(cfg.elem,1),1);
% cfg.vessel = vessel;
% cfg.radius = 2*ones(size(cfg.elem,1),1);
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

%% Compare vessel MMC with and without connection point

% 1, 18, 22
n = 22;

figure,
load /drives/neza2/users/yaoshen/NEU/Research/mmc/without_connection.mat
subplot(121),imagesc(log10(squeeze(flux.data(1:50,1:60,n))),[-4 8]);
colorbar
title(['Without connection, slice=' num2str(n)]),colormap jet

load /drives/neza2/users/yaoshen/NEU/Research/mmc/with_connection.mat
subplot(122),imagesc(log10(squeeze(flux.data(1:50,1:60,n))),[-4 8]);
colorbar
title(['With connection, slice=' num2str(n)]),colormap jet

%% Compare original MMC and vessel MMC

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