% run iMMC on 99th vessel

clear
% close all

addpath(genpath('/space/neza/2/users/yaoshen/NEU/Research/Redbird/tensorlab/'))
addpath(genpath('/space/neza/2/users/yaoshen/NEU/Research/iso2mesh/'))

%% Split element with more than 2 vessel edges

load /drives/neza2/users/yaoshen/NEU/Research/mmc/mmc/vessel99

centroid = elemcentroid(node_tet,elem_tet);
elemindex = find(sum(localedge,2)>2);

n2e = [1 2; 1 3; 1 4; 2 3; 2 4; 3 4];
faceorder = [1 2 3; 1 2 4; 1 3 4; 2 3 4];
nindex = size(node_tet,1);
node_tet = [node_tet; zeros(length(elemindex),3)];
newelem = zeros(4,5,length(elemindex));
newvessel = zeros(4,6,length(elemindex));
for i=1:length(elemindex)
    cen = centroid(elemindex(i),:);
    nindex = nindex+1;
    node_tet(nindex,:) = cen;
    et = elem_tet(elemindex(i),:);
    clear vesselnode
    vesselnode = {'null'};
    for ie=1:6
        if localedge(elemindex(i),ie)
            vesselnode{end+1} = num2str(et(n2e(ie,:)));
            vesselnode{end+1} = num2str(fliplr(et(n2e(ie,:))));
        end
    end
    for j=1:4
         newelemt = [et(faceorder(j,:)) nindex et(5)];
         newelem(j,:,i) = newelemt;
        for ie=1:6
            edgenew = num2str(newelemt(n2e(ie,:)));
            if ismember(edgenew,vesselnode)
                newvessel(j,ie,i) = 1;
            end
        end
    end
end

new_elem = elem_tet(1:elemindex(1)-1,:);
new_localedge = localedge(1:elemindex(1)-1,:);
for i=1:length(elemindex)-1
    post = elem_tet(elemindex(i)+1:elemindex(i+1)-1,:);
    new_elem = [new_elem; newelem(:,:,i); post];
    
    post2 = localedge(elemindex(i)+1:elemindex(i+1)-1,:);
    new_localedge = [new_localedge; newvessel(:,:,i); post2];
end
post = elem_tet(elemindex(end)+1:end,:);
new_elem = [new_elem; newelem(:,:,end); post];

post2 = localedge(elemindex(end)+1:end,:);
new_localedge = [new_localedge; newvessel(:,:,end); post2];

%% Check if splitting is successful

% figure,plotmesh(node_tet,new_elem,'x>50')
find(sum(new_localedge,2)>2)

%% Plot vessel from mesh (check the label is correct)

plotvessel = 0;
if plotvessel==1
    n2e = [1 2; 1 3; 1 4; 2 3; 2 4; 3 4];
    stack{1} = 'haha';
    figure,
    for i=1:size(new_localedge,1)
        if sum(new_localedge(i,:))
            plotmesh(node_tet,new_elem(i,:),'facealpha',0,'facecolor','w','edgealpha',0.1);
        end
        for j=1:6
            if new_localedge(i,j)
                edget = [new_elem(i,n2e(j,1)) new_elem(i,n2e(j,2))];
                edge_str = num2str(sort(edget));
                if ~ismember(edge_str,stack)
                    stack{end+1} = edge_str;

    %                 edget
                    xe = [node_tet(edget(1),1) node_tet(edget(2),1)];
                    ye = [node_tet(edget(1),2) node_tet(edget(2),2)];
                    ze = [node_tet(edget(1),3) node_tet(edget(2),3)];
                    plot3(xe,ye,ze,'r','LineWidth',max([noder(edget(1)) noder(edget(2))]))
    %                 plot3(xe,ye,ze,'r')
                    hold on
                end
            end
        end
    end
    axis equal
end

%%
% new_elem = elem_tet;
% new_localedge = localedge;
vessel = 7*ones(size(new_localedge,1),2);
radius = zeros(size(new_localedge,1),2);

for i=1:size(vessel,1)
    index = find(new_localedge(i,:)==1);
    if ~isempty(index)
        vessel(i,1) = index(1);
        edget1 = [new_elem(i,n2e(index(1),1)) new_elem(i,n2e(index(1),2))];
        radius(i,1) = max([noder(edget1(1)) noder(edget1(2))]);
        if length(index)>1
            vessel(i,2) = index(2);
            edget2 = [new_elem(i,n2e(index(2),1)) new_elem(i,n2e(index(2),2))];
            radius(i,2) = max([noder(edget2(1)) noder(edget2(2))]);
        end
    end
end

elem = [new_elem(:,1:4) vessel-1 6*ones(size(vessel)) radius zeros(size(vessel))];

%%
% elem = [elem_tet(:,1:4) 6*ones(size(elem_tet,1),4) zeros(size(elem_tet,1),4)];
% elem = [new_elem(:,1:4) vessel-1 6*ones(size(vessel)) 0.01*ones(size(vessel)) zeros(size(vessel))];
%% run mmc

clear cfg
cfg.nphoton=1e6;
cfg.node = [node_tet zeros(size(node_tet,1),1)];
cfg.elem = elem;
cfg.implicit = 1;
cfg.elemprop=ones(size(cfg.elem,1),1);
% cfg.vessel = vessel;
% cfg.radius = 2*ones(size(cfg.elem,1),1);
cfg.srcpos=[50.5 50.5 0];
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

% figure,plotmesh(cfg.node,cfg.elem(9,[1 2 3]),'facealpha',0.2),hold on
% plotmesh(cfg.node,cfg.elem(9,[1 2 4]),'facealpha',0.2)
% plotmesh(cfg.node,cfg.elem(9,[1 3 4]),'facealpha',0.2)
% plotmesh(cfg.node,cfg.elem(9,[2 3 4]),'facealpha',0.2)
% xlabel('x'),ylabel('y'),zlabel('z')
% 
% plotmesh(cfg.node,cfg.elem(115,[1 2 3]),'facealpha',0.2),hold on
% plotmesh(cfg.node,cfg.elem(115,[1 2 4]),'facealpha',0.2)
% plotmesh(cfg.node,cfg.elem(115,[1 3 4]),'facealpha',0.2)
% plotmesh(cfg.node,cfg.elem(115,[2 3 4]),'facealpha',0.2)

flux=mmclab(cfg);
% save flux_vessel_blood flux

figure,slice3(log10(flux.data))
caxis([-5 8])
colormap jet

%% Element centroid

function centroid = elemcentroid(node,elem)

centroid = zeros(size(elem,1),3);
for i=1:size(elem,1)
    n1 = node(elem(i,1),:);
    n2 = node(elem(i,2),:);
    n3 = node(elem(i,3),:);
    n4 = node(elem(i,4),:);
    n = [n1; n2; n3; n4];
    centroid(i,:) = mean(n);
end

end