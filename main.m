clear all
close all 
clc

for time_ind = 0:100000
    filename = ['ExportCellProp_' num2str(time_ind) '.txt'];
    display(filename);
    temp = exist(filename);
    while(temp==0)
       pause(2);
       temp = exist(filename);
    end
    epi_nodes = load(['ExportCellProp_' num2str(time_ind) '.txt']);
    while(size(epi_nodes,1)==0)
       pause(2);
       epi_nodes = load(['ExportCellProp_' num2str(time_ind) '.txt']);
    end
    while(epi_nodes(end,1)~=123456789)
       pause(2);
       epi_nodes = load(['ExportCellProp_' num2str(time_ind) '.txt']);
    end
    main_signaling(time_ind); 
end
   



