function [newcell] = sort_counterclock(newcell,centroid)
angle = atan((newcell(:,2)-centroid(2))./(newcell(:,1)-centroid(1))).*(newcell(:,1)>=centroid(1))+...
    (pi+atan((newcell(:,2)-centroid(2))./(newcell(:,1)-centroid(1)))).*(newcell(:,1)<centroid(1));
angle(angle<0) = angle(angle<0) + 2*pi;
[val,ind] = sort(angle,'ascend');  
newcell = newcell(ind,:);