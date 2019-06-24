% Main test function for verifying the projection method for calculating
% the intersection between the ray and the vessel cylinder
clear
% close all

% define the cylinder
E0 = [0 0 -20];
E1 = [0 0 10];
radius = 5;

% generate random P0 and P1 in a cube
N = 100;
range = [-10 10];
P0 = (range(2)-range(1))*rand(N,3)+range(1);
P1 = (range(2)-range(1))*rand(N,3)+range(1);

tic
count = 1;
for i=1:N
    % get the hitting point
    H = line_cylinder_intersect(P0(i,:),P1(i,:),E0,E1,radius);
    if sum(isinf(H))>0
        continue;
    end
    HP(count,:) = [H i];
    count = count+1;
end
t = toc

tic
count = 1;
for i=1:N
    % get the hitting point
    H2 = line_cylinder_intersect2(P0(i,:),P1(i,:),E0,E1,radius);
    if sum(isinf(H2))>0
        continue;
    end
    HP2(count,:) = [H2 i];
    count = count+1;
end
t2 = toc

figure,plot3(HP(:,1),HP(:,2),HP(:,3),'b.','MarkerSize',16,'LineWidth',3)
hold on
for j=1:count-1
    dp = P1(HP(j,4),:)-P0(HP(j,4),:);
%     plot3([P0(HP(j,4),1); P1(HP(j,4),1)],[P0(HP(j,4),2); P1(HP(j,4),2)],[P0(HP(j,4),3); P1(HP(j,4),3)]);
        quiver3(P0(HP(j,4),1),P0(HP(j,4),2),P0(HP(j,4),3),dp(1),dp(2),dp(3));
end
% plot3(HP2(:,1),HP2(:,2),HP2(:,3),'r.')
xlim([-10 10]),ylim([-10 10]),zlim([-10 10])
xlabel('x'),ylabel('y'),zlabel('z')
axis equal

plot3(HP2(:,1),HP2(:,2),HP2(:,3),'r.','MarkerSize',16,'LineWidth',3)
% hold on
% for j=1:count-1
%    plot3([P0(HP2(j,4),1); P1(HP2(j,4),1)],[P0(HP2(j,4),2); P1(HP2(j,4),2)],[P0(HP2(j,4),3); P1(HP2(j,4),3)]) 
% end
% % plot3(HP2(:,1),HP2(:,2),HP2(:,3),'r.')
% xlim([-10 10]),ylim([-10 10]),zlim([-10 10])
% axis equal
% xlabel('x'),ylabel('y'),zlabel('z')