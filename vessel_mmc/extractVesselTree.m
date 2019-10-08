% Use DFS to obtain tree structure from the centerline data

clear all
% close all

addpath(genpath('/space/neza/2/users/yaoshen/NEU/Research/Redbird/tensorlab/'))
addpath('/space/neza/2/users/yaoshen/NEU/Research/Optimal_Speration/Tutorial_MMC-MCX/nii')
addpath('/space/neza/2/users/yaoshen/NEU/Research/Optimal_Speration/Tutorial_MMC-MCX/nfri_functions')
addpath('/space/neza/2/users/yaoshen/NEU/Research/Optimal_Speration/Tutorial_MMC-MCX/metch')

%% Read vessel center line and radius

vesselr = load_nii('/space/wazu/3/users/yaoshen/backup/iMMC/Vessel-data/radius/99.nii.gz');
r = vesselr.img;
% figure,slice3(r(100:200,100:200,100:200))

vesselc = load_nii('/space/wazu/3/users/yaoshen/backup/iMMC/Vessel-data/centerline/99.nii.gz');
center = vesselc.img;
% figure,slice3(center(100:200,100:200,100:200))

vesselb = load_nii('/space/wazu/3/users/yaoshen/backup/iMMC/Vessel-data/bif/99.nii.gz');
bif = vesselb.img;
% figure,slice3(bif(100:200,100:200,100:200))

load bifcenter_99.mat % bifcenter bradii

%% Use DSF to obtain vessel tree structure

cc = 1;
for i=-1:1
    for j=-1:1
        for k=-1:1
            if i==0 && j==0 && k==0
                continue;
            end
            dir(cc,:) = [i j k];
            cc = cc+1;
        end
    end
end

r_trim = r; % (100:200,100:200,100:200);
bi_trim = bif; %(100:200,100:200,100:200);
[bx,by,bz] = size(r_trim);

global rvol bvol bindex bifctr segcount segment index count STACK
bvol = bi_trim;
bindex = bradii;
bifctr = bifcenter;
rvol = r_trim;
segcount = 0;
segment = 1000; % save the node every 'segment'
index = 1;
count = 1;
STACK(1).x=0; STACK(1).y=0; STACK(1).z=0; STACK(1).radius=0; STACK(1).index=0; STACK(1).parent=NaN;% STACK(1).children=zeros(30,1);

for zz=10:size(rvol,3)
    slice = squeeze(rvol(:,:,zz));
    while ~isempty(find(slice~=0,1))
        [rr,cc] = find(slice~=0);
        root.x=rr(1); root.y=cc(1); root.z=zz;
        root.radius = slice(rr(1),cc(1));
        root.index = index;
        index = index+1;
        root.parent = NaN;
        STACK(count) = root;
        count = count+1;
        
        rvol(root.x,root.y,root.z) = 0;
        dsf_tracevessel(root,dir);
        slice = squeeze(rvol(:,:,zz));
    end
end

%% Plot vessel segments

for i=1:length(STACK)
    xx(i) = STACK(i).x;
    yy(i) = STACK(i).y;
    zz(i) = STACK(i).z;
    rr(i) = STACK(i).radius;
    parent(i) = STACK(i).parent;
end
figure,scatter3(xx,yy,zz,rr,'r','filled')
xlabel('x'),ylabel('y'),zlabel('z')
title(['segment=' num2str(segment)])
axis equal

%% Simplify tree

nc = 1;
prev = [xx(1) yy(1) zz(1)];
noder(nc) = rr(10);
node(nc,:) = prev;
node_str{nc} = num2str(prev);
nc = nc+1;
for i=2:length(STACK)
    curr = [xx(i) yy(i) zz(i)];
    curr_str = num2str(curr);    
    if ~ismember(curr_str,node_str)
        noder(nc) = rr(i);
        node(nc,:) = curr;
        node_str{nc} = curr_str;
        nc = nc+1;
    end
end

% new bifurcation
b2n = zeros(length(STACK),1);
for i=1:length(STACK)
    xid=find(node(:,1)==xx(i)); yid=find(node(:,2)==yy(i)); zid=find(node(:,3)==zz(i));
    nindex = intersect(intersect(xid,yid),zid);
    if ~isempty(nindex)
        b2n(i) = nindex;
    end
end

% new parent
parentt = parent;
b2nt = [b2n; NaN];
parentt(isnan(parentt)) = length(b2nt);
parent_new = b2nt(parentt);

%% extract edge

ne = 1;
edge(ne,:) = [b2n(2) parent_new(2)];
edge_str{ne} = num2str(edge(ne,:));
for i=2:length(b2n)
    if isnan(parent_new(i))
        continue;
    end
    if b2n(i)==parent_new(i)
        continue;
    end
    curedge_str = num2str([b2n(i) parent_new(i)]);
    if ~ismember(curedge_str,edge_str)
        ne = ne+1;
        edge(ne,:) = [b2n(i) parent_new(i)];
        edge_str{ne} = num2str(edge(ne,:));
    end
end
fprintf('haha')

%%
figure,
for i=1:size(edge,1)
    i
    xe = [node(edge(i,1),1) node(edge(i,2),1)];
    ye = [node(edge(i,1),2) node(edge(i,2),2)];
    ze = [node(edge(i,1),3) node(edge(i,2),3)];
    plot3(xe,ye,ze,'r','LineWidth',noder(edge(i,1)))
    hold on
    axis equal
end

%% Squeeze stack

% end_flag = 0;
% ns = 1;
% bcount = 1;
% while ~end_flag
%     num_start = b2n(ns);
%     b2n_squeeze(bcount) = num_start;
%     while b2n(ns)==num_start
%         if ns+1>length(b2n)
%             end_flag = 1;
%             break;
%         end
%         ns = ns+1;
%     end
%     bcount = bcount+1;
% end


%% Get bifurcation in vessel segments

% pu = unique(parent);
% inan = find(isnan(pu)==1,1);
% pu = pu(1:inan-1);
% N = histc(parent,pu);
% bi_perent = pu(N>=2);

%%
function dsf_tracevessel(current,dir)
    global rvol bindex bifctr segcount segment index count STACK
    [bx,by,bz] = size(rvol);
    
    for nb=1:size(dir,1)
        neighbor = current;
        neighbor.x = neighbor.x+dir(nb,1);
        neighbor.y = neighbor.y+dir(nb,2);
        neighbor.z = neighbor.z+dir(nb,3);

        if neighbor.x<=0 || neighbor.x>bx || neighbor.y<=0 || neighbor.y>by...
                || neighbor.z<=0 || neighbor.z>bz
           continue; 
        end
        if rvol(neighbor.x,neighbor.y,neighbor.z)==0
            continue;
        end
        
        neighbor.radius = rvol(neighbor.x,neighbor.y,neighbor.z);
        segcount = segcount+1;
        
        bifidx = bindex(neighbor.x,neighbor.y,neighbor.z);
        if segcount>segment
            neighbor.index = index;
            index = index+1;
            neighbor.parent = current.index;
            segcount = 0;
            STACK(count) = neighbor;
            if bifidx>0
                STACK(count).x=bifctr(bifidx,1); STACK(count).y=bifctr(bifidx,2); STACK(count).z=bifctr(bifidx,3);
            end
            count = count+1;
        elseif bindex(neighbor.x,neighbor.y,neighbor.z)>0
            neighbor.index = index;
            index = index+1;
            neighbor.parent = current.index;
            STACK(count) = neighbor;
            STACK(count).x=bifctr(bifidx,1); STACK(count).y=bifctr(bifidx,2); STACK(count).z=bifctr(bifidx,3);
            count = count+1;
        end

        rvol(neighbor.x,neighbor.y,neighbor.z) = 0;
        dsf_tracevessel(neighbor,dir);
        
    end
    
    return;
end