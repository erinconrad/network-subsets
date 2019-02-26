function [overlap,dist] = compare_elecs(group_1,group_2,locs,doPlot)

% This function takes 2 groups of electrodes and compares how close they
% are together

% get electrode locations
locs_1 = locs(group_1,:);
locs_2 = locs(group_2,:);

% get centroids
centroid_1 = mean(locs_1,1);
centroid_2 = mean(locs_2,1);

% get the distance between centroids
dist = sqrt(sum((centroid_1-centroid_2).^2));

% get the percentage of electrodes in GROUP 1 that are also in group 2
num_overlap = sum(ismember(group_1,group_2));
overlap = num_overlap/length(group_1);

if doPlot == 1
    figure
    scatter3(locs(:,1),locs(:,2),locs(:,3),100,'k');
    hold on
    scatter3(locs_1(:,1),locs_1(:,2),locs_1(:,3),50,'r','filled');
    hold on
    scatter3(locs_2(:,1),locs_2(:,2),locs_2(:,3),20,'b','filled');
    scatter3(centroid_1(1),centroid_1(2),centroid_1(3),200,'r');
    scatter3(centroid_2(1),centroid_2(2),centroid_2(3),200,'b');
    fprintf('The number overlapping is %d and the distance is %1.1f.\n',...
        num_overlap,dist);
    pause
    close(gcf);
    
end

end