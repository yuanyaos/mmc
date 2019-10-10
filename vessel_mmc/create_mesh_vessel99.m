clear
close all

load /drives/neza2/users/yaoshen/NEU/Research/mmc/mmc/vessel_mmc/vessel99
thre = 300;
[node,edge,noder] = reducevesselsize(node,edge,noder,thre,thre,thre);
figure,plotvessel(edge,node,noder)

[nbox,ebox] = meshgrid6(0:1,0:1,0:1);
fbox = volface(ebox);
nbox = nbox*(thre+20)-10;
plotbox(nbox,ebox)

% [nbox,fbox,ebox] = meshabox([0 0 0],[300 300 300],100000,100000);
% hold on
% plotmesh(nbox,fbox,'facealpha',0.1,'edgealpha',0.2,'facecolor','c')

noffset = size(nbox,1);
edge = edge+noffset;
noder = [zeros(1,noffset) noder];
node = [nbox; node];
[nn,nd] = size(node);

fedge = [edge edge(:,2)];
face = [fbox; fedge];

sessionid = 'vessel99';
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

%% Read and plot vessel mesh

clear

elem = importdata('/space/neza/2/users/yaoshen/NEU/Research/mmc/mmc/vessel_mmc/tetgen/vessel99/vessel99.1.ele');
face = importdata('/space/neza/2/users/yaoshen/NEU/Research/mmc/mmc/vessel_mmc/tetgen/vessel99/vessel99.1.face');
node = importdata('/space/neza/2/users/yaoshen/NEU/Research/mmc/mmc/vessel_mmc/tetgen/vessel99/vessel99.1.node');

elem = elem(2:end,2:end);
elem = elem+1;
face = face(2:end,2:end);
face = face+1;
node = node(2:end,2:end);

elem = [elem ones(size(elem,1),1)];

figure,plotmesh(node,elem,'x>100','facealpha',1)


%% Add function

function [node_new,edge_new,noder_new] = reducevesselsize(node,edge,noder,thresholdx,thresholdy,thresholdz)
    count = 1;
    for i=1:size(edge,1)
        if node(edge(i,1),1)>thresholdx || node(edge(i,2),1)>thresholdx || ...
                node(edge(i,1),2)>thresholdy || node(edge(i,2),2)>thresholdy || ...
                node(edge(i,1),3)>thresholdz || node(edge(i,2),3)>thresholdz
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
    plot3(xe,ye,ze,'r','LineWidth',noder(edge(i,1)))
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