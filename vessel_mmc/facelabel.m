function [faceelemid, elem_edge] = facelabel(elem,node,face,ploton)
% elem: mesh element
% node: mesh node
% face: face of the surface
% ploton: 1/0->plot/not plot
%
% faceelemid: indices of elements that contain the surface
% elem_edge: save the reference element indices for the elements that
% contain all surface edges but not surface faces

if nargin<4
    ploton = 0;
end

% find elements that contain the surface face
faceelemid = [];
for i=1:size(face,1)
    [row1,~] = find(elem==face(i,1));
    [row2,~] = find(elem==face(i,2));
    [row3,~] = find(elem==face(i,3));
    
    faceelemid = [faceelemid; intersect(intersect(row1,row2),row3)];
end
faceelemid = unique(faceelemid);

% extract all edges of the surface
edge = [];
for i=1:size(face,1)
	e1(1) = face(i,1); % edge 1
    e1(2) = face(i,2);
    
    e2(1) = face(i,1); % edge 2
    e2(2) = face(i,3);
    
    e3(1) = face(i,2); % edge 3
    e3(2) = face(i,3);
    
    edge = [edge; e1; e2; e3];
end
edge = unique(edge,'rows');

% find elements that contain the edges
elem_edge = zeros(size(elem,1),1);
for i=1:size(edge,1)
    [row1,~] = find(elem==edge(i,1));
    [row2,~] = find(elem==edge(i,2));
    
    edfcelem = intersect(row1,row2);    % elements that contain the edge
    edgeelem = setdiff(edfcelem,faceelemid);    % elements that contain the edge but is not contained in face element
    elemref = intersect(edfcelem,faceelemid);   % elements that contain the edge but is also contained in face element
    
    elem_edge(edgeelem) = elemref(1);   % record the reference element for the edge
end


%% Plot
if ploton
    index = (elem_edge~=0);
    figure,plotmesh(node,elem(faceelemid,:),'facealpha',0.1,'edgealpha',0.5,'facecolor','b')
    hold on    
    plotmesh(node,elem(index,:),'facealpha',0.1,'edgealpha',0.5,'facecolor','r')
    
    figure,
%     for i=1:length(elem_edge)
%         if(elem_edge(i)==0)
%             continue;
%         end
%         plotmesh(node,elem(i,:),'facealpha',0.1,'edgealpha',0.5,'facecolor','r')
%         hold on
%         scatter3(node(elem(i,:),1),node(elem(i,:),2),node(elem(i,:),3))
%         plotmesh(node,elem(elem_edge(i),:),'facealpha',0.1,'edgealpha',0.5,'facecolor','g')
%     end
    plotmesh(node,elem(index,:),'facealpha',0.1,'edgealpha',0.5,'facecolor','r')
    hold on
    extra_index = nonzeros(elem_edge);
    plotmesh(node,elem(extra_index,:),'facealpha',0.1,'edgealpha',0.5,'facecolor','g')
end

end