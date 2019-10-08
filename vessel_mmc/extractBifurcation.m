% Use DFS to obtain tree structure from the centerline data

clear all
% close all

addpath(genpath('/space/neza/2/users/yaoshen/NEU/Research/Redbird/tensorlab/'))
addpath('/space/neza/2/users/yaoshen/NEU/Research/Optimal_Speration/Tutorial_MMC-MCX/nii')
addpath('/space/neza/2/users/yaoshen/NEU/Research/Optimal_Speration/Tutorial_MMC-MCX/nfri_functions')
addpath('/space/neza/2/users/yaoshen/NEU/Research/Optimal_Speration/Tutorial_MMC-MCX/metch')

%% Read vessel center line and radius

vesselb = load_nii('/space/wazu/3/users/yaoshen/backup/iMMC/Vessel-data/bif/99.nii.gz');
bif = vesselb.img;
% figure,slice3(bif(100:200,100:200,100:200))

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

bif_trim = bif; %(100:200,100:200,100:200);
[bx,by,bz] = size(bif_trim);

global bvol bradii index
bvol = bif_trim;
bradii = zeros(size(bvol));
segcount = 0;
segment = 10; % save the node every 'segment'
index = 1;

for zz=1:size(bvol,3)
    zz
    slice = squeeze(bvol(:,:,zz));
    while ~isempty(find(slice~=0,1))
        [rr,cc] = find(slice~=0);
        root.x=rr(1); root.y=cc(1); root.z=zz;
        bradii(root.x,root.y,root.z) = index;        
        bvol(root.x,root.y,root.z) = 0;
        dsf_tracevessel(root,dir);
        
        slice = squeeze(bvol(:,:,zz));
        index = index+1;
    end
    
end

%%

bcenter = zeros(size(bradii));
for i=1:index-1
    i
    [r,c,v] = ind2sub(size(bradii),find(bradii==i));
    bifcenter(i,:) = [mean(r) mean(c) mean(v)];
    bifradius(i,:) = bif_trim(r(1),c(1),v(1));
    bcenter(round(bifcenter(i,1)),round(bifcenter(i,2)),round(bifcenter(i,3))) = i;
end
figure,slice3(bradii)

%%
function dsf_tracevessel(current,dir)
    global bvol bradii index
    [bx,by,bz] = size(bvol);
    
    for nb=1:size(dir,1)
        neighbor = current;
        neighbor.x = neighbor.x+dir(nb,1);
        neighbor.y = neighbor.y+dir(nb,2);
        neighbor.z = neighbor.z+dir(nb,3);

        if neighbor.x<=0 || neighbor.x>bx || neighbor.y<=0 || neighbor.y>by...
                || neighbor.z<=0 || neighbor.z>bz
           continue; 
        end
        if bvol(neighbor.x,neighbor.y,neighbor.z)==0
            continue;
        end

        bradii(neighbor.x,neighbor.y,neighbor.z) = index;
        bvol(neighbor.x,neighbor.y,neighbor.z) = 0;
        dsf_tracevessel(neighbor,dir);
        
    end
    
    return;
end