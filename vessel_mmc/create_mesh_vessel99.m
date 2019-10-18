clear
close all

%% Create node and face for vessel

load vessel99
% ====== give the range of bounding box ======
thre = [100 200];


[node,edge,noder] = reducevesselsize(node,edge,noder,thre,thre,thre);
figure,plotvessel(edge,node,noder)

[nbox,ebox] = meshgrid6(0:1,0:1,0:1);
fbox = volface(ebox);
nbox = nbox*(thre(2)-thre(1)+10)+thre(1)-5;
plotbox(nbox,ebox)

% [nbox,fbox,ebox] = meshabox([-5 -5 -5],[205 205 205],10000,10000);
% hold on
% plotmesh(nbox,fbox,'facealpha',0.1,'edgealpha',0.2,'facecolor','c')

noffset = size(nbox,1);
edge = edge+noffset;
noder = [zeros(1,noffset) noder];
node = [nbox; node];
[nn,nd] = size(node);

% fedge = [edge edge(:,2)];
fedge = [edge(:,1) edge];
face = [fbox; fedge];

%% Write .poly file for tetgen

% ====== save file ======
savePoly = 1;

sizestr = '100';
sessionid = ['vessel99_' sizestr];

if savePoly    
    fid=fopen([sessionid,'.poly'],'wt');
    % === node ===
    fprintf(fid,'#node list\n');
    fprintf(fid,'%d %d %d %d\n',nn,nd,0,0);
    fprintf(fid,'%d %6.16f %6.16f %6.16f\n',[0:nn-1; node(:,1)'; node(:,2)'; node(:,3)']);

    % === face ===
    [ne,~] = size(face);
    face = face-1;
    fprintf(fid,'#face list\n');
    fprintf(fid,'%d %d\n',ne,1);
    for i=1:ne
        fprintf(fid,'%d %d\n',1,0);
        fprintf(fid,'%d %d %d %d\n',3,face(i,1),face(i,2),face(i,3));
    end

    % === hole ===
    fprintf(fid,'#hole list\n');
    fprintf(fid,'%d\n',0);

    fclose(fid);
end

%% Read and plot vessel mesh

% clear
path = ['/space/neza/2/users/yaoshen/NEU/Research/mmc/mmc/vessel_mmc/tetgen/' sessionid '/' sessionid '.1.ele'];
elem_tet = importdata(path);
path = ['/space/neza/2/users/yaoshen/NEU/Research/mmc/mmc/vessel_mmc/tetgen/' sessionid '/' sessionid '.1.face'];
face_tet = importdata(path);
path = ['/space/neza/2/users/yaoshen/NEU/Research/mmc/mmc/vessel_mmc/tetgen/' sessionid '/' sessionid '.1.node'];
node_tet = importdata(path);

elem_tet = elem_tet(2:end,2:end);
elem_tet = elem_tet+1;
face_tet = face_tet(2:end,2:end);
face_tet = face_tet+1;
node_tet = node_tet(2:end,2:end);

elem_tet = [elem_tet ones(size(elem_tet,1),1)];

figure,plotmesh(node_tet,elem_tet,'x>50','facealpha',1)

%% Label vessel in mesh

% compare node before/after tetgen
ctet = 1;
nid = [];
for i=1:min([size(node,1) size(node_tet,1)])
    if ~isequal(node(i,:),node_tet(i,:))
        nid(ctet) = i;
        ctet = ctet+1;
    end
end

if isempty(nid)     % if the order of node does not change
    localedge = labelvessel(elem_tet,edge);
else
    error('The order of nodes before and after tetgen does not match')
end

%% Plot vessel from mesh (check the label is correct)

n2e = [1 2; 1 3; 1 4; 2 3; 2 4; 3 4];
stack{1} = 'haha';
figure,
for i=1:size(localedge,1)
    if sum(localedge(i,:))
        plotmesh(node_tet,elem_tet(i,:),'facealpha',0,'facecolor','w','edgealpha',0.1);
    end
    for j=1:6
        if localedge(i,j)
            edget = [elem_tet(i,n2e(j,1)) elem_tet(i,n2e(j,2))];
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


edge_missed = checkedges(localedge,edge,elem_tet);

for ie=2:length(edge_missed)
    em = str2num(edge_missed{ie});
    xe = [node_tet(em(1),1) node_tet(em(2),1)];
    ye = [node_tet(em(1),2) node_tet(em(2),2)];
    ze = [node_tet(em(1),3) node_tet(em(2),3)];
    plot3(xe,ye,ze,'r','LineWidth',noder(em(1)))
    hold on
end




%% Functions

function [node_new,edge_new,noder_new] = reducevesselsize(node,edge,noder,thresholdx,thresholdy,thresholdz)
    count = 1;
    for i=1:size(edge,1)
        if node(edge(i,1),1)<thresholdx(1) || node(edge(i,2),1)<thresholdx(1) || node(edge(i,1),1)>thresholdx(2) || node(edge(i,2),1)>thresholdx(2) || ...
           node(edge(i,1),2)<thresholdy(1) || node(edge(i,2),2)<thresholdy(1) || node(edge(i,1),2)>thresholdy(2) || node(edge(i,2),2)>thresholdy(2) || ...
           node(edge(i,1),3)<thresholdz(1) || node(edge(i,2),3)<thresholdz(1) || node(edge(i,1),3)>thresholdz(2) || node(edge(i,2),3)>thresholdz(2)
           continue;
        end
        edge_new(count,:) = edge(i,:);
        count = count+1;
    end
    
    noid = unique(edge_new(:));
    nnmap = zeros(size(node,1),1);
    for j=1:length(noid)
        nnmap(noid(j)) = j;
        node_new(j,:) = node(noid(j),:);
        noder_new(j) = noder(noid(j));
    end
    edge_new = nnmap(edge_new);
end

%% Plot vessel

function plotvessel(edge,node,noder)

for i=1:size(edge,1)
%     i
    xe = [node(edge(i,1),1) node(edge(i,2),1)];
    ye = [node(edge(i,1),2) node(edge(i,2),2)];
    ze = [node(edge(i,1),3) node(edge(i,2),3)];
    plot3(xe,ye,ze,'r','LineWidth',max([noder(edge(i,1)) noder(edge(i,2))]))
%     plot3(xe,ye,ze,'r')
    hold on
end
axis equal

end

%% Plot bounding box

function plotbox(node,elem)

face = volface(elem);
plotmesh(node,face,'facealpha',0.05,'facecolor','c')
% for i=1:size(elem,1)
%     plotmesh(node,[elem(i,1) elem(i,2) elem(i,3)],'facealpha',0.05,'facecolor','c')
%     hold on
%     plotmesh(node,[elem(i,1) elem(i,2) elem(i,4)],'facealpha',0.05,'facecolor','c')
%     plotmesh(node,[elem(i,1) elem(i,3) elem(i,4)],'facealpha',0.05,'facecolor','c')
%     plotmesh(node,[elem(i,2) elem(i,3) elem(i,4)],'facealpha',0.05,'facecolor','c')
% end
axis equal

end

%% Label local vessel edge for every element

function localedge = labelvessel(elem,edge)
% node: node coordinates (nn x 3)
% elem: node indices for every element (ne x 4)
% edge: node indice for vessel edge (nedge x 2)
%
% localedge: local vessel edge index for every element (ne x 6)

elem = elem(:,1:4);
% elem_ori = elem;

localedge = zeros(size(elem,1),6);

stack_edge = {'hahaha'};
for i=1:size(edge,1)
    stack_edge{end+1} = num2str(edge(i,:));
    stack_edge{end+1} = num2str([edge(i,2) edge(i,1)]);
end

% elem = sort(elem,2);
% edge = sort(edge,2);
order = [1 2; 1 3; 1 4; 2 3; 2 4; 3 4];
for i=1:size(elem,1)
    for j=1:6   % 6 edges in each tetrahedron
        edgelocal = [elem(i,order(j,1)) elem(i,order(j,2))];
        edgelocal_str = num2str(edgelocal);
        if ismember(edgelocal_str,stack_edge)
            localedge(i,j) = 1;
        end
%         le = edge(:,1)==edgelocal(1);
%         ri = edge(:,2)==edgelocal(2);
%         if sum(le & ri)>0   % if edgelocal is vessel, then obtain its local edge index 
%             for k=1:6
%                 edge_ori = [elem_ori(i,order(k,1)) elem_ori(i,order(k,2))];
%                 if isempty(setdiff(edge_ori,edgelocal))
%                     localedge(i,k) = 1;
%                 end
%             end
%         end
    end
end

end

%% Check edges
% Check if original edges match edges after tetgen

function stack_edge = checkedges(localedge,edge,elem_tet)
% edge: edges that should be contained in new element after tetgen
% elem_tet: new element after tetgen

n2e = [1 2; 1 3; 1 4; 2 3; 2 4; 3 4];
stack_edge{1} = 'haha';
stack_localedge{1} = 'haha';

% for i=1:size(edge,1)
%     stack_edge{end+1} = num2str(edge(i,:));
% %     stack_edge{end+1} = num2str([edge(i,2) edge(i,1)]);
% end

for i=1:size(localedge,1)
    for j=1:6
        if localedge(i,j)
            edge_str = num2str([elem_tet(i,n2e(j,1)) elem_tet(i,n2e(j,2))]);
            edge_str2 = num2str([elem_tet(i,n2e(j,2)) elem_tet(i,n2e(j,1))]);
            if ~ismember(edge_str,stack_localedge)
                stack_localedge{end+1} = edge_str;
            end
            if ~ismember(edge_str2,stack_localedge)
                stack_localedge{end+1} = edge_str2;
            end
        end
    end
end

for i=1:size(edge,1)
    edge_str_ori = num2str(edge(i,:));
    if ~ismember(edge_str_ori,stack_localedge)
        stack_edge{end+1} = edge_str_ori;
    end
end

end