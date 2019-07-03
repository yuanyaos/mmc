function [vessel,vesseln,pathelem,nodeelem] = vessellabel(elem,node,pse,radius,ploton)
% pse(1): node index of starting point
% pse(2): node index of ending point

% vessel: nex1
% vesseln: nnx1
% noderadius: nnx1

addpath('/drives/neza2/users/yaoshen/NEU/Research/mmc/vessel_mmc/dijkstra')

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

ps = find(eleuniq==pse(1));
pe = find(eleuniq==pse(2));
[costs,paths] = dijkstra(A,C,ps,pe);

pindex = eleuniq(paths);
pnode = node(pindex,:);

%% get vessel index
edgenodelocalindex = [];
nodelocalindex = [];
pathelem = [];
nodeelem = [];
[row1,col1] = find(elem==pindex(1));
for i=1:length(pindex)-1
    [row2,col2] = find(elem==pindex(i+1));
    [elemt,ia,ib] = intersect(row1,row2);
    [nodet,in] = setdiff(row2,elemt);
    edgenodelocalindex = [edgenodelocalindex; [col1(ia) col2(ib)]-1];
    nodelocalindex = [nodelocalindex; col2(in)-1];
    pathelem = [pathelem; elemt];
    nodeelem = [nodeelem; nodet];
    row1 = row2;
    col1 = col2;
end

% map from local node index to local edge index
n2e = {'0  1', '0  2', '0  3', '1  2', '1  3', '2  3'};
vessel = 6*ones(size(elem,1),1);
vesseln = zeros(size(node,1),1);
[pathelem,ie,~] = unique(pathelem);
edgenodelocalindex = edgenodelocalindex(ie,:);
edgenodelocalindex = sort(edgenodelocalindex,2);
[nodeelem,ie,~] = unique(nodeelem);
nodelocalindex = nodelocalindex(ie,:);
for i=1:size(edgenodelocalindex,1)
    vt(i) = find(strcmp(n2e,num2str(edgenodelocalindex(i,:))))-1;
    vessel(pathelem(i)) = vt(i);
end
% vessel = vessel+1;

for i=1:length(nodeelem)
    vesseln(elem(nodeelem(i),nodelocalindex(i)+1)) = radius;
end

%% Plot
if ploton
    figure,plotmesh(node,elem,'facealpha',0.1,'edgealpha',0.5,'facecolor','w')
    hold on
    plot3(pnode(:,1),pnode(:,2),pnode(:,3),'LineWidth',2)
    plotmesh(node,elem(pathelem,:),'facealpha',0.1,'edgealpha',0.5,'facecolor','r')
    plotmesh(node,elem(nodeelem,:),'facealpha',0.1,'edgealpha',0.5,'facecolor','b')
    idx = find(vesseln~=0);
    plot3(node(idx,1),node(idx,2),node(idx,3),'og','LineWidth',3);
    xlabel('x'),ylabel('y'),zlabel('z')
end

end