function plotvessel(elem,node)
% Plot edge vessel and node vessel

n2e = {'0  1', '0  2', '0  3', '1  2', '1  3', '2  3'};

vessel = elem(:,5:6);
vesselr = elem(:,7:8);
elem = elem(:,1:4);

noder = node(:,4);
node = node(:,1:3);

figure,plotmesh(node,elem,'facealpha',0.1,'edgealpha',0.5,'facecolor','w')
hold on

for i=1:size(vessel,1)
    for j=1:size(vessel,2)
        if vessel(i,j)==6
            continue;
        end
        v = str2num(n2e{vessel(i,j)+1});
        node1 = node(elem(i,v(1)+1),:);
        node2 = node(elem(i,v(2)+1),:);
        node12 = [node1; node2];
        plot3(node12(:,1),node12(:,2),node12(:,3),'r','LineWidth',2*vesselr(i,j))
    end
end

for i=1:size(node,1)
    if noder(i)
        scatter3(node(i,1),node(i,2),node(i,3),50*noder(i),'b','filled')
    end
end
xlabel('x'),ylabel('y'),zlabel('z')

end

