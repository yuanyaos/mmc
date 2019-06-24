clear
addpath('/space/neza/2/users/yaoshen/NEU/Research/iso2mesh/')
addpath('/drives/neza2/users/yaoshen/NEU/Research/mmc/vessel_mmc/dijkstra')

%%

[node,face,elem]=meshabox([0 0 0],[60 70 80],1000,1000);
% plotmesh(node,elem);

elem = elem(:,1:4);
eleuniq = unique(elem(:));
A = zeros(length(eleuniq));
for i=1:length(eleuniq)
    [row,~] = find(elem==eleuniq(i));
    t = unique(elem(row,:));
    A(i,t) = 1;
end
A = A-eye(size(A));

C = zeros(size(A));
for i=1:length(eleuniq)
    for j=i+1:length(eleuniq)
        C(i,j) = norm(node(eleuniq(i),:)-node(eleuniq(j),:));
    end
end
C = C+C';

%%
close all

% define start and end points
pse = round(rand(2,1)*length(eleuniq));
% [costs,paths] = dijkstra(A,C,pse(1),pse(2));
[costs,paths] = dijkstra(A,C,1,181);

pindex = eleuniq(paths);
pnode = node(pindex,:);

figure,plotmesh(node,elem,'facealpha',0.1,'edgealpha',0.5,'facecolor','w')
hold on
plot3(pnode(:,1),pnode(:,2),pnode(:,3),'-o','LineWidth',2)
plot3(pnode(1,1),pnode(1,2),pnode(1,3),'-og','LineWidth',2)
plot3(pnode(end,1),pnode(end,2),pnode(end,3),'-or','LineWidth',2)
xlabel('x'),ylabel('y'),zlabel('z')

%% get vessel index

nodelocalindex = [];
elempath = [];
[row1,col1] = find(elem==pindex(1));
for i=1:length(pindex)-1   
    [row2,col2] = find(elem==pindex(i+1));
    [elemt,ia,ib] = intersect(row1,row2);
    nodelocalindex = [nodelocalindex; [col1(ia) col2(ib)]-1];
    elempath = [elempath; elemt];
    row1 = row2;
    col1 = col2;
end

% map from local node index to local edge index
n2e = {'0  1', '0  2', '0  3', '1  2', '1  3', '2  3'};
vessel = zeros(size(elem,1),1);
[elempath,ie,~] = unique(elempath);
nodelocalindex = nodelocalindex(ie,:);
for i=1:size(nodelocalindex,1)
    nodelocalindex = sort(nodelocalindex,2);
    vt(i) = find(strcmp(n2e,num2str(nodelocalindex(i,:))))-1;
    vessel(elempath(i)) = vt(i);
end

elemvessel = elem(elempath,:);
plotmesh(node,elemvessel,'facealpha',0.1,'edgealpha',0.5,'facecolor','r')

%% run vessel mmc

addpath(genpath('/drives/neza2/users/yaoshen/NEU/Research/mmc/mmc'))
addpath(genpath('/space/neza/2/users/yaoshen/NEU/Research/Redbird/tensorlab/'))

clear cfg
cfg.nphoton=1e7;
cfg.node = node;
cfg.elem = elem;
cfg.elemprop=ones(size(cfg.elem,1),1);
cfg.vessel = vessel;
cfg.radius = 2*ones(size(cfg.elem,1),1);
cfg.srcpos=[2 2 0];
cfg.srcdir=[0 0 1];
cfg.prop=[0 0 1 1;0.005 1 0 1.37;0.5 1 0 1.37];
cfg.tstart=0;
cfg.tend=5e-9;
cfg.tstep=5e-9;
cfg.isreflect=1;
cfg.debuglevel='TP';
cfg.issaveref=1;  % in addition to volumetric fluence, also save surface diffuse reflectance
cfg.method = 'grid';

% run the simulation
flux=mmclab(cfg);

figure,slice3(log10(flux.data))