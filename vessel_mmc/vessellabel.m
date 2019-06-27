function [vessel] = vessellabel(elem,node,pse)
% pse(1): node index of starting point
% pse(2): node index of ending point

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
vessel = 6*ones(size(elem,1),1);
[elempath,ie,~] = unique(elempath);
nodelocalindex = nodelocalindex(ie,:);
for i=1:size(nodelocalindex,1)
    nodelocalindex = sort(nodelocalindex,2);
    vt(i) = find(strcmp(n2e,num2str(nodelocalindex(i,:))))-1;
    vessel(elempath(i)) = vt(i);
end

end