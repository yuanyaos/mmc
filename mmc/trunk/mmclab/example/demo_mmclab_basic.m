%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MMCLAB - Mesh-based Monte Carlo for MATLAB/Octave by Qianqina Fang
%
% In this example, we show the most basic usage of MMCLAB.
%
% This file is part of Mesh-based Monte Carlo (MMC) URL:http://mcx.sf.net/mmc
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% prepare simulation input

cfg.nphoton=1e6;
[cfg.node face cfg.elem]=meshabox([0 0 0],[60 60 30],6);
cfg.elemprop=ones(size(cfg.elem,1),1);
cfg.srcpos=[30 30 0];
cfg.srcdir=[0 0 1];
cfg.prop=[0 0 1 1;0.005 1 0 1.37];
cfg.tstart=0;
cfg.tend=5e-9;
cfg.tstep=5e-9;
cfg.debuglevel='TP';

% run the simulation

flux=mmclab(cfg);

% plotting the result

% if you have the SVN version of iso2mesh, use the next line to plot:
% qmeshcut(cfg.elem(:,1:4),cfg.node(:,1:3),log10(abs(flux.data(:))),'y=30','linestyle','none');

[cutpos,cutvalue,facedata]=qmeshcut(cfg.elem(:,1:4),cfg.node(:,1:3),log10(abs(flux.data(:))),'y=30');
hcut=patch('Faces',facedata,'Vertices',cutpos,'FaceVertexCData',cutvalue,'facecolor','interp','linestyle','none');
view([0 1 0]);
set(gca,'dataaspectratio',[1 1 1]);
colorbar;
