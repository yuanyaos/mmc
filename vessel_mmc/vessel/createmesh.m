sessionid='vessel';

addpath(genpath('/space/neza/2/users/yaoshen/NEU/Research/Redbird/tensorlab/'))
addpath(genpath('/space/neza/2/users/yaoshen/NEU/Research/iso2mesh/'))
addpath('/drives/neza2/users/yaoshen/NEU/Research/mmc/mmc/vessel_mmc/dijkstra')
addpath('/drives/neza2/users/yaoshen/NEU/Research/mmc/mmc/matlab/');

[node,elem] = meshgrid6(0:1,0:1,0:1);
node = node*60;
node = [node [1 0 0 0 0 0 0 1]'];
pse(1) = 1;
pse(2) = 8;
% [vessel] = vessellabel(elem,node,pse);
[vessel,vesseln,pathelem,nodeelem] = vessellabel(elem,node,pse,2,0);
elem = [elem ones(size(vessel)) vessel 6*ones(size(vessel)) ones(size(vessel)) ones(size(vessel))];

% elem=sortrows(elem);
% elem(:,1:4)=meshreorient(node,elem(:,1:4));

srcpos=[30.5 30.5 0];
savemmcmesh(sessionid,node,elem,[]);
eid=tsearchn(node(:,1:3),elem(:,1:4),srcpos)

%% write optical property

prop = [0 0 1 1;0.005 1 0 1.37;0.5 1 0 1.37];
fid=fopen(['prop_',sessionid,'.dat'],'wt');
fprintf(fid,'%d\t%d\n',1,size(prop,1)-1);
fprintf(fid,'%d\t%d\t%d\t%d\t%d\n', [1:size(prop,1)-1;prop(2:end,:)']);
fclose(fid);

%% plot

addpath(genpath('/space/neza/2/users/yaoshen/NEU/Research/Redbird/tensorlab/'))

filename = 'vessel.dat';
fid=fopen(filename,'rt');
hd=fscanf(fid,'%f');
fclose(fid);

fluence = hd(2:2:end);
fluence = reshape(fluence,[61 61 61]);

figure,slice3(log10(fluence))
caxis([-5 8])
colormap jet