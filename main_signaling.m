function main_signaling(index_couple)

global thres_lateral thres_edge thres_nodes
thres_lateral = 0.64; % need to adjust carefully to avoid triangles at junctions
thres_edge = 0.005;
thres_nodes = 0.64;

% load cell configuration from mechanical submodel
epi_nodes = load(['ExportCellProp_' num2str(index_couple) '.txt']);
epi_nodes(end,:) = [];

total_cell = max(epi_nodes(:,1))+1;
cell_nodes = epi_nodes(total_cell+1:end,:);
cell_nodes = [cell_nodes(:,2:3), cell_nodes(:,1)];
for i = 1:total_cell
    eval(['centroid_' num2str(i-1) '=epi_nodes(i,2:3);']);
end

dist_nodes = pdist2(cell_nodes(:,1:2),cell_nodes(:,1:2));
% indicator of detection of common edge between neighboring cells
common_edge_id = zeros(total_cell,total_cell);

for i = 0:total_cell-1
    row_i = find(cell_nodes(:,3)==i);
    dist_i = dist_nodes(row_i,:);
    [row,col] = find(dist_i<thres_lateral);
    neighbor_i = unique(cell_nodes(col,3));
    neighbor_i = neighbor_i(neighbor_i~=i);
    %construct common edge and indicator for cell-cell contact
    cella = cell_nodes(cell_nodes(:,3)==i,1:2);
    for j = 1:length(neighbor_i)
        if common_edge_id(i+1,neighbor_i(j)+1)==0
            cellb = cell_nodes(cell_nodes(:,3)==neighbor_i(j),1:2);
            eval(['[intf_' num2str(i) '_' num2str(neighbor_i(j)) ', index_' num2str(i) '_' num2str(neighbor_i(j)) ', index_' num2str(neighbor_i(j)) '_' num2str(i) '] = detect_common_edge(cella,cellb);']);
            common_edge_id(i+1,neighbor_i(j)+1)=1;
            common_edge_id(neighbor_i(j)+1,i+1)=1;
        end
    end
end

%figure(1); hold on;
%plot(cell_nodes(:,1),cell_nodes(:,2),'r.');

% plot(intf_0_1(:,1),intf_0_1(:,2),'o-');
% plot(intf_0_4(:,1),intf_0_4(:,2),'o-');
% plot(intf_0_5(:,1),intf_0_5(:,2),'o-');
% plot(centroid_0(:,1),centroid_0(:,2),'o');
% plot(centroid_1(:,1),centroid_1(:,2),'o');
% plot(centroid_4(:,1),centroid_4(:,2),'o');
% plot(centroid_5(:,1),centroid_5(:,2),'o'); 

% detect vertices on each cell one by one
for i = 0:total_cell-1
    id = find(cell_nodes(:,3)==i);
    newcell = cell_nodes(id,1:2);
    neighbor_i = find(common_edge_id(i+1,:));
    neighbor_i = neighbor_i - 1;
    delete_edge = [];
    add_edge = [];
    for j = 1:length(neighbor_i)
        eval(['delete_edge = [delete_edge; index_' num2str(i) '_' num2str(neighbor_i(j)) '];']);
        if i<neighbor_i(j)
            eval(['add_edge = [add_edge; intf_' num2str(i) '_' num2str(neighbor_i(j)) '];']);
        else
            eval(['add_edge = [add_edge; intf_' num2str(neighbor_i(j)) '_' num2str(i) '];']);
        end
    end 
    newcell(delete_edge,:) = [];
    newcell = [newcell; add_edge];
            
    % sort the nodes counter-clockwisely
    eval(['[newcell] = sort_counterclock(newcell,centroid_' num2str(i) ');']);
%     plot(newcell(:,1),newcell(:,2),'*'); pause

    eval(['[vt_' num2str(i) '] = detect_vertices(newcell);']);
end

% plot(vt_1(:,1),vt_1(:,2),'*','Color',[0 0 0]); pause

%choose the union set of vertices on common edge
for i = 0:total_cell-1
    neighbor_i = find(common_edge_id(i+1,:));
    neighbor_i = neighbor_i - 1;
    for j = 1:length(neighbor_i)
        if i<neighbor_i(j)
            eval(['[C_keep,ia_keep,ib_keep] = intersect(vt_' num2str(neighbor_i(j)) ',intf_' num2str(i) '_' num2str(neighbor_i(j)) ',''rows'');']);
        else
            eval(['[C_keep,ia_keep,ib_keep] = intersect(vt_' num2str(neighbor_i(j)) ',intf_' num2str(neighbor_i(j)) '_' num2str(i) ',''rows'');']);
        end
        eval(['vt_' num2str(i) '= [vt_' num2str(i) ';C_keep];']);
    end
    eval(['[vt_' num2str(i) '] = sort_counterclock(vt_' num2str(i) ',centroid_' num2str(i) ');']);
    eval(['[vt_' num2str(i) '] = unique(vt_' num2str(i) ',''rows'');']);
    
end

% for i = 0:total_cell-1
%     eval(['plot(vt_' num2str(i) '(:,1),vt_' num2str(i) '(:,2),''o'');']); 
% end


%replace points too close by the middle point on common edge
for i = 0:total_cell-1
    neighbor_i = find(common_edge_id(i+1,:));
    neighbor_i = neighbor_i - 1;
    for j = 1:length(neighbor_i)
        if neighbor_i(j)>i
            eval(['[vt_contact,ia,ib] = intersect(vt_' num2str(i) ',vt_' num2str(neighbor_i(j)) ',''rows'');']);
            eval(['[vt_contact] = sort_counterclock(vt_contact,centroid_' num2str(i) ');']); 
        
            eval(['[vt_contact_' num2str(i) '_' num2str(neighbor_i(j)) '] = combine_nodes(vt_contact);']);
             
            eval(['vt_' num2str(i) '(ia,:) = [];']); 
            eval(['vt_' num2str(i) ' = [vt_' num2str(i) ';vt_contact];']);
            eval(['vt_' num2str(neighbor_i(j)) '(ib,:) = [];']);
            eval(['vt_' num2str(neighbor_i(j)) '= [vt_' num2str(neighbor_i(j)) ';vt_contact];']);
            eval(['[vt_' num2str(i) '] = sort_counterclock(vt_' num2str(i) ',centroid_' num2str(i) ');']);
            eval(['[vt_' num2str(neighbor_i(j)) '] = sort_counterclock(vt_' num2str(neighbor_i(j)) ',centroid_' num2str(neighbor_i(j)) ');']);
        end
    end        
end

% replace junction points among multiple(>=3) contacting cells by centroids
for i = 0:total_cell-1
    vt_num = eval(['size(vt_' num2str(i) ',1)']);
    eval(['vt_dist = sqrt((vt_' num2str(i) '(:,1)-[vt_' num2str(i) '(2:end,1); vt_' num2str(i) '(1,1)]).^2+(vt_' num2str(i) '(:,2)-[vt_' num2str(i) '(2:end,2); vt_' num2str(i) '(1,2)]).^2);']);
    [val,ind] = min(vt_dist);
    while val<thres_nodes
        if ind < vt_num
            eval(['pt1 = vt_' num2str(i) '(ind,:);']);
            eval(['pt2 = vt_' num2str(i) '(ind+1,:);']);
            eval(['mid_pt = ( vt_' num2str(i) '(ind,:)+vt_' num2str(i) '(ind+1,:) )/2;']);
            eval(['vt_' num2str(i) '(ind,:) = mid_pt;']);
            eval(['vt_' num2str(i) '(ind+1,:) = [];']);
        else
            eval(['pt1 = vt_' num2str(i) '(end,:);']);
            eval(['pt2 = vt_' num2str(i) '(1,:);']);
            eval(['mid_pt = ( vt_' num2str(i) '(end,:)+vt_' num2str(i) '(1,:) )/2;']);
            eval(['vt_' num2str(i) '(end,:) = mid_pt;']);
            eval(['vt_' num2str(i) '(1,:) = [];']);
        end
        neighbor_i = find(common_edge_id(i+1,:));
        neighbor_i = neighbor_i - 1;
        for j = 1:length(neighbor_i)
            eval(['[id1,ind1] = ismember(pt1,vt_' num2str(neighbor_i(j)) ',''rows'');']);
            if id1 == 1
                eval(['vt_' num2str(neighbor_i(j)) '(ind1,:) = mid_pt;']);
            end
            eval(['[id2,ind2] = ismember(pt2,vt_' num2str(neighbor_i(j)) ',''rows'');']);
            if id2 == 1
                eval(['vt_' num2str(neighbor_i(j)) '(ind2,:) = mid_pt;']);
            end
        end
        vt_num = vt_num -1;
        eval(['vt_dist = sqrt((vt_' num2str(i) '(:,1)-[vt_' num2str(i) '(2:end,1); vt_' num2str(i) '(1,1)]).^2+(vt_' num2str(i) '(:,2)-[vt_' num2str(i) '(2:end,2); vt_' num2str(i) '(1,2)]).^2);']);
        [val,ind] = min(vt_dist);
    end
end
        
% visualize the triangular mesh
%for i = 0:total_cell-1
%    eval(['plot(vt_' num2str(i) '(:,1),vt_' num2str(i) '(:,2),''*'');']); 
%end
%saveas(figure(1),'vertices.fig'); 

for i = 0:total_cell-1
    eval(['temp = vt_' num2str(i) ';']);
%    plot([temp(:,1); temp(1,1)],[temp(:,2); temp(1,2)],'b')
%    for j = 1:size(temp,1)
%        eval(['plot([centroid_' num2str(i) '(1) temp(' num2str(j) ',1)],[centroid_' num2str(i) '(2) temp(' num2str(j) ',2)],''b'');']);
%    end
end
%saveas(figure(1),'triangle_mesh.fig'); 


% solve u_t = D(u_xx+u_yy)-du+source on the triangular mesh by
% approximating diffusion by passive transport and imposing source function
% at cells close to 25
mesh_size = 0;
vt_all = [];
for i = 0:total_cell-1
    eval(['num_tri_' num2str(i) ' = size(vt_' num2str(i) ',1);']);
    eval(['vt_all = [vt_all; vt_' num2str(i) '];'])
    eval(['mesh_size = mesh_size + num_tri_' num2str(i) ';']);
end
u = zeros(mesh_size,1);
% construct a matrix to denote neighboring triangles
NI_mat = [];
% construct a matrix to store contact length
A_mat = [];
% intracellular neighbors
for i = 0:total_cell-1
    eval(['temp = ones(num_tri_' num2str(i) ',1);']);
    eval(['NI_i = spdiags(temp,1,num_tri_' num2str(i) ',num_tri_' num2str(i) ')+spdiags(temp,-1,num_tri_' num2str(i) ',num_tri_' num2str(i) ');']);
    eval(['NI_i(1,num_tri_' num2str(i) ') = 1;']);
    eval(['NI_i(num_tri_' num2str(i) ',1) = 1;']);
    NI_mat = blkdiag(NI_mat,NI_i);
    
    eval(['edge_vec = pdist2(centroid_' num2str(i) ',vt_' num2str(i) ');']);
    edge_vec = [edge_vec(2:end), edge_vec(1)];
    edge_vec = edge_vec';
    eval(['A_' num2str(i) ' = spdiags(edge_vec,1,num_tri_' num2str(i) ',num_tri_' num2str(i) ') + spdiags(edge_vec,-1,num_tri_' num2str(i) ',num_tri_' num2str(i) ');']);
    eval(['A_' num2str(i) '(1,end) = edge_vec(end);']);
    eval(['A_' num2str(i) '(end,1) = edge_vec(end);']);
    eval(['A_mat = blkdiag(A_mat,A_' num2str(i) ');']); 
end
% intercellular neighbors
for i = 0:total_cell-1
    eval(['num_tri = num_tri_' num2str(i) ';']);
    for j = 1:num_tri
        if j == num_tri
            jnext = 1;
        else
            jnext = j+1;
        end
        eval(['[ind1] = find(sum(abs(vt_all-repmat(vt_' num2str(i) '(j,:),mesh_size,1)),2)==0);']);
        eval(['[ind2] = find(sum(abs(vt_all-repmat(vt_' num2str(i) '(jnext,:),mesh_size,1)),2)==0);']);
        if length(ind1)>1 && length(ind2)>1
            eval(['Ledge = pdist(vt_' num2str(i) '([j jnext],:));']);
            ind1_sort = sort(ind1);
            ind2_sort = sort(ind2); 
            if abs(ind1_sort(1)-ind2_sort(1)) ==1
                tri_1 = min(ind1_sort(1),ind2_sort(1));
            else
                tri_1 = max(ind1_sort(1),ind2_sort(1));
            end
            if abs(ind1_sort(2)-ind2_sort(2)) ==1
                tri_2 = min(ind1_sort(2),ind2_sort(2));
            else
                tri_2 = max(ind1_sort(2),ind2_sort(2));
            end
            NI_mat(tri_1,tri_2) = 1;
            NI_mat(tri_2,tri_1) = 1;
            A_mat(tri_1,tri_2) = Ledge;
            A_mat(tri_2,tri_1) = Ledge;
        end
    end
end

% centroids of each triangle
mesh_centroids = zeros(mesh_size,2);
num_tri_pre = 0;
for i = 0:total_cell-1
    eval(['num_tri = num_tri_' num2str(i) ';']);
    for j = 1:num_tri
        if j == num_tri
            jnext = 1;
        else
            jnext = j+1;
        end
        eval(['mesh_centroids(num_tri_pre+j,:) =  sum([vt_' num2str(i) '([j jnext],:);centroid_' num2str(i) '])/3;']);
    end
    num_tri_pre = num_tri_pre + num_tri; 
end
L_mat = pdist2(mesh_centroids,mesh_centroids);
L_mat = L_mat + 10000*eye(mesh_size); %avoid singularities

% source cells: located within 12% of the total tissue size around the
% midline
tissue_centroid = zeros(1,2);
for i = 1:total_cell
    eval(['tissue_centroid = tissue_centroid + centroid_' num2str(i-1) ';']);
end
tissue_centroid = tissue_centroid/total_cell;
tissue_r = sqrt( sum((tissue_centroid-centroid_0).^2) );
for i = 2:total_cell
    eval(['temp = sqrt( sum((tissue_centroid-centroid_' num2str(i-1) ').^2) );']);
    if temp > tissue_r
        tissue_r = temp;
    end
end
    
[ind] = find(abs(cell_nodes(:,1)-tissue_centroid(1))<0.12*tissue_r);
source_cellid = unique(cell_nodes(ind,3));
source_cellid = sort(source_cellid,'ascend');
u_source = [];
for i = 0:total_cell-1
    temp = ismember(i,source_cellid);
    if temp == 1
        eval(['u_source = [u_source; ones(num_tri_' num2str(i) ',1)];']);
    else
        eval(['u_source = [u_source; zeros(num_tri_' num2str(i) ',1)];']);
    end
end


Dpp_mat = zeros(mesh_size,4);
dt = 0.002;

for iter = 1:1000
    [frhs] = Dpp_signaling(Dpp_mat,A_mat,L_mat,NI_mat,u_source);
    temp = Dpp_mat + dt * frhs;
    
    [frhs] = Dpp_signaling(temp,A_mat,L_mat,NI_mat,u_source);
    Dpp_mat = Dpp_mat/2 + temp/2 + dt/2 * frhs;
end

Dpp = Dpp_mat(:,1);
Dpp_cell = zeros(total_cell,1);
Tkv = Dpp_mat(:,2);
Tkv_cell = zeros(total_cell,1);
DT = Dpp_mat(:,3);
DT_cell = zeros(total_cell,1);
pMad = Dpp_mat(:,4);
pMad_cell = zeros(total_cell,1);
%figure(2); 
%subplot(2,2,1); hold on;
num_tri_pre=0;
for i = 0:total_cell-1
    eval(['num_tri = num_tri_' num2str(i) ';']);
    for j = 1:num_tri
        if j == num_tri
            jnext = 1;
        else
             jnext = j+1;
         end
%         eval(['fill([vt_' num2str(i) '(j,1) vt_' num2str(i) '(jnext,1) centroid_' num2str(i) '(1)],[vt_' num2str(i) '(j,2) vt_' num2str(i) '(jnext,2) centroid_' num2str(i) '(2)],Dpp(num_tri_pre+j));']);
     end
     Dpp_cell(i+1) = sum(Dpp(num_tri_pre+1:num_tri_pre+num_tri))/num_tri;
     num_tri_pre = num_tri_pre + num_tri;
 end
% colorbar;
% subplot(2,2,2); hold on;
 num_tri_pre=0;
 for i = 0:total_cell-1
     eval(['num_tri = num_tri_' num2str(i) ';']);
     for j = 1:num_tri
         if j == num_tri
             jnext = 1;
         else
             jnext = j+1;
         end
%         eval(['fill([vt_' num2str(i) '(j,1) vt_' num2str(i) '(jnext,1) centroid_' num2str(i) '(1)],[vt_' num2str(i) '(j,2) vt_' num2str(i) '(jnext,2) centroid_' num2str(i) '(2)],Tkv(num_tri_pre+j));']);
     end
     Tkv_cell(i+1) = sum(Tkv(num_tri_pre+1:num_tri_pre+num_tri))/num_tri;
     num_tri_pre = num_tri_pre + num_tri;
 end
% colorbar;
% subplot(2,2,3); hold on;
 num_tri_pre=0;
 for i = 0:total_cell-1
     eval(['num_tri = num_tri_' num2str(i) ';']);
     for j = 1:num_tri
         if j == num_tri
             jnext = 1;
         else
             jnext = j+1;
         end
%         eval(['fill([vt_' num2str(i) '(j,1) vt_' num2str(i) '(jnext,1) centroid_' num2str(i) '(1)],[vt_' num2str(i) '(j,2) vt_' num2str(i) '(jnext,2) centroid_' num2str(i) '(2)],DT(num_tri_pre+j));']);
     end
     DT_cell(i+1) = sum(DT(num_tri_pre+1:num_tri_pre+num_tri))/num_tri;
     num_tri_pre = num_tri_pre + num_tri;
 end
% colorbar;
% subplot(2,2,4); hold on;
 num_tri_pre=0;
 for i = 0:total_cell-1
     eval(['num_tri = num_tri_' num2str(i) ';']);
     for j = 1:num_tri
         if j == num_tri
             jnext = 1;
         else
             jnext = j+1;
         end
%         eval(['fill([vt_' num2str(i) '(j,1) vt_' num2str(i) '(jnext,1) centroid_' num2str(i) '(1)],[vt_' num2str(i) '(j,2) vt_' num2str(i) '(jnext,2) centroid_' num2str(i) '(2)],pMad(num_tri_pre+j));']);
     end
     pMad_cell(i+1) = sum(pMad(num_tri_pre+1:num_tri_pre+num_tri))/num_tri;
     num_tri_pre = num_tri_pre + num_tri;
 end
% colorbar;
% saveas(figure(2),'Dpp_signaling.fig');

filename = ['Dpp_cell_T' num2str(index_couple) '.txt'];
%fileID = fopen(filename,'w');
fileID = fopen('tmp.txt','w+');
fprintf(fileID,'%.4f\n',Dpp_cell);
fclose(fileID);
movefile ('tmp.txt',filename) ; 




