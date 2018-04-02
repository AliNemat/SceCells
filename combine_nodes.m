function [compact_nodes] = combine_nodes(nodes)
global thres_nodes
compact_nodes = nodes;
points = compact_nodes(1:end-1,:);
points_next = compact_nodes(2:end,:);
dist_neighbor = sqrt((points(:,1)-points_next(:,1)).^2+(points(:,2)-points_next(:,2)).^2);
while min(dist_neighbor)<thres_nodes
    [val,ind] = min(dist_neighbor);
    insert_point = (points(ind,:)+points_next(ind,:))/2;
    compact_nodes(ind,:) = insert_point;
    compact_nodes(ind+1,:) = [];
    points = compact_nodes(1:end-1,:);
    points_next = compact_nodes(2:end,:);
    dist_neighbor = sqrt((points(:,1)-points_next(:,1)).^2+(points(:,2)-points_next(:,2)).^2);
end