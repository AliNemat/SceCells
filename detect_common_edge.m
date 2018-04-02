function [intf,id_cell12,id_cell21] = detect_common_edge(cell1,cell2)
global thres_lateral
D = pdist2(cell1,cell2);
[D1,I1] = min(D,[],2);
[D2,I2] = min(D,[],1);
[id_cell12]=find(D1<thres_lateral);
[id_cell21]=find(D2<thres_lateral);
id_cell21 = id_cell21';
intf = [(cell1(id_cell12,1)+cell2(I1(id_cell12),1))/2 (cell1(id_cell12,2)+cell2(I1(id_cell12),2))/2];

